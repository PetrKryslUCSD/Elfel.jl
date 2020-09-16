"""
    klein_gordon_wavelet

The Klein–Gordon equation is often written in natural units:
```latex
psi_tt - psi_xx + m^2 psi = 0
```

Here it is solved on a finite interval, with essential boundary condition in the
form of a bell function. Only the real part of psi is solved for.
"""
module klein_gordon_wavelet

using LinearAlgebra
using MeshCore.Exports
using MeshSteward
using MeshSteward.Exports
using Elfel.Exports
using PlotlyJS

L = 20.0
m2 = 5.0 # coefficient m^2
N = 100;# number of subdivisions along the length of the domain

# Boundary and initial conditions
Fbc(x) = 0.0
Fic(x) = exp(-1.5*x^2)
Vic(x) = 0.0
tend = 5.0
dt = 0.01

function genmesh()
    ir = L2block(L, N)
    MeshSteward.transform(ir, x -> x-L/2)
    return attach!(Mesh(), ir)
end

function assembleM(Psih)
    function integrate!(am, geom, elit, qpit)
        nedof = ndofsperel(elit)
        me = LocalMatrixAssembler(nedof, nedof, 0.0)
        for el in elit
            init!(me, eldofs(el), eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                N = bfun(qp)
                for j in 1:nedof
                    for i in 1:nedof
                        me[i, j] += (N[i] * N[j]) * (JxW)
                    end
                end
            end
            assemble!(am, me)
        end
        return am
    end

    elit = FEIterator(Psih)
    qpit = QPIterator(Psih, (kind = :Gauss, order = 2))
    geom = geometry(Psih.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Psih), ndofs(Psih))

    integrate!(am, geom, elit, qpit)

    return finish!(am)
end

function assembleK(Psih)
    function integrate!(am, geom, elit, qpit)
        nedof = ndofsperel(elit)
        me = LocalMatrixAssembler(nedof, nedof, 0.0)
        for el in elit
            init!(me, eldofs(el), eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                x = location(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                N = bfun(qp)
                for j in 1:nedof
                    for i in 1:nedof
                        me[i, j] += (gradN[i][1] * gradN[j][1]) * (JxW)
                    end
                end
            end
            assemble!(am, me)
        end
        return am
    end

    elit = FEIterator(Psih)
    qpit = QPIterator(Psih, (kind = :Gauss, order = 2))
    geom = geometry(Psih.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Psih), ndofs(Psih))

    integrate!(am, geom, elit, qpit)

    return finish!(am)
end

function assembleD(Psih, m2)
    function integrate!(am, geom, elit, qpit, m)
        nedof = ndofsperel(elit)
        me = LocalMatrixAssembler(nedof, nedof, 0.0)
        for el in elit
            init!(me, eldofs(el), eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                x = location(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                N = bfun(qp)
                for j in 1:nedof
                    for i in 1:nedof
                        me[i, j] += (N[i] * N[j]) * (m2 * JxW)
                    end
                end
            end
            assemble!(am, me)
        end
        return am
    end

    elit = FEIterator(Psih)
    qpit = QPIterator(Psih, (kind = :Gauss, order = 2))
    geom = geometry(Psih.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Psih), ndofs(Psih))

    integrate!(am, geom, elit, qpit, m2)

    return finish!(am)
end

function run()
    mesh = genmesh()
    Psih = FESpace(Float64, mesh, FEH1_L2())
    ir = baseincrel(mesh)
    bdvl = connectedv(boundary(mesh));
    locs = geometry(mesh)
    for i in bdvl
        setebc!(Psih, 0, i, 1, Fbc(locs[i]...))
    end
    numberdofs!(Psih)
    @show nu = nunknowns(Psih)
    xs = [locs[i][1] for i in 1:length(locs)]
    M = assembleM(Psih)
    K = assembleK(Psih)
    D = assembleD(Psih, m2)
    Psi0 = fill(0.0, ndofs(Psih))
    V0 = fill(0.0, ndofs(Psih))
    for i in 1:length(locs)
        n = dofnum(Psih, 0, i, 1)
        Psi0[n] = Fic(locs[i][1])
        V0[n] = Vic(locs[i][1])
    end
    @show Psi0
    for i in bdvl
        n = dofnum(Psih, 0, i, 1)
        V0[n] = 0.0
    end
    V1 = fill(0.0, ndofs(Psih))
    L1 = fill(0.0, ndofs(Psih))
    L0 = fill(0.0, ndofs(Psih))
    # Plots
    layout = Layout(;width=700, height=700, xaxis = attr(title="x", range=[-L/2, L/2]), yaxis = attr(title = "Psi", range=[-1.0, 1.0]))
    sigdig = n -> round(n*1000)/1000
    function updateplot(pl, t, xs, Psis)
        curv = scatter(;x=xs, y=Psis, mode="lines", name = "Sol", line_color = "rgb(155, 15, 15)")
        plots = cat(curv; dims = 1)
        pl.plot.layout["title"] = "t = $(sigdig(t))"
        react!(pl, plots, pl.plot.layout)
        sleep(0.12)
    end
    function solve!(V1, M, K, D, L0, L1, Psi0, V0, dt, nu)
        # This assumes zero EBC
        R = (L1 + L0)/2 + (1/dt)*(M*V0) - K*(Psi0 + (dt/4)*V0) - D*(Psi0 + (dt/4)*V0)
        V1[1:nu] = ((1/dt)*M + (dt/4)*(K + D))[1:nu, 1:nu] \ R[1:nu]
        return V1
    end
    pl = 0
    step = 0
    t = 0.0
    while t < tend
        if mod(step, 10) == 0
            Psis = [Psi0[dofnum(Psih, 0, i, 1)] for i in 1:nshapes(ir.right)]
            if step == 0
                curv = scatter(;x=xs, y=Psis, mode="lines")
                plots = cat(curv; dims = 1)
                pl = plot(plots, layout)
                display(pl)
            end
            updateplot(pl, t, xs, Psis)
        end
        solve!(V1, M, K, D, L0, L1, Psi0, V0, dt, nu)
        Psi1 = Psi0 + dt/2*(V1 + V0)
        t = t + dt
        step = step + 1
        copyto!(V0, V1)
        copyto!(Psi0, Psi1)
    end
    # scattersysvec!(Psih, Psi0)
    # makeattribute(Psih, "Psi", 1)
    # vtkwrite("klein_gordon_wavelet-Psi", baseincrel(mesh), [(name = "Psi",)])
end

end

using .klein_gordon_wavelet
klein_gordon_wavelet.run()
