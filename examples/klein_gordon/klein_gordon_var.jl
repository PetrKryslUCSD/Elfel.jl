"""
    klein_gordon_var

The Klein–Gordon equation with a variable coefficient:
```latex
psi_tt - c^2 psi_xx + g psi = 0
```
where ``g = g(t, x)``.

Here it is solved on a finite interval, with zero essential boundary conditions at the ends of the interval.

Only the real part of psi is solved for.
"""
module klein_gordon_var

using LinearAlgebra
using MeshCore.Exports
using MeshSteward.Exports
using Elfel.Exports
using PlotlyJS

L = 2.0
c = 10.0 # coefficient
N = 150;# number of subdivisions along the length of the domain
Fbc(x) = 0.0

# G(t, x) = 800*sin(5*pi*t)*(x/L)^3*(L-x)/L;# 

# Harmonic vibration
g(t, x) = 0.0; 
Fic(x) = 0.0
Vic(x) = 5*sin(pi*x/L) - 7*sin(3*pi*x/L)
tend = 80.0
dt = 0.05

# Uniform velocity
g(t, x) = 0.0;
Fic(x) = 0.0
Vic(x) = 7.0
tend = 2.0
dt = 0.002

# Uniform velocity, nonzero constant G
g(t, x) = -100;
Fic(x) = 0.0
Vic(x) = 7.0
tend = 2.0
dt = 0.002

# # Uniform velocity, nonzero G variable in x
g(t, x) = -1000.0*(x/L)^3;
Fic(x) = 0.0
Vic(x) = -7.0
tend = 2.0
dt = 0.001

# Uniform velocity, nonzero G variable in t
g(t, x) = +200*sin(1.32*pi*t) # 
Fic(x) = 0.0
Vic(x) = 7*sin(pi*x/L)
tend = 18.0
dt = 0.005

function genmesh()
    return attach!(Mesh(), L2block(L, N))
end

function assembleM(Fh)
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

    elit = FEIterator(Fh)
    qpit = QPIterator(Fh, (kind = :Gauss, order = 2))
    geom = geometry(Fh.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Fh), ndofs(Fh))

    integrate!(am, geom, elit, qpit)

    return finish!(am)
end

function assembleK(Fh, c)
    function integrate!(am, geom, elit, qpit, c)
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
                        me[i, j] += (gradN[i][1] * gradN[j][1]) * (c^2 * JxW)
                    end
                end
            end
            assemble!(am, me)
        end
        return am
    end

    elit = FEIterator(Fh)
    qpit = QPIterator(Fh, (kind = :Gauss, order = 2))
    geom = geometry(Fh.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Fh), ndofs(Fh))

    integrate!(am, geom, elit, qpit, c)

    return finish!(am)
end

function assembleD(Fh, g, t)
    function integrate!(am, geom, elit, qpit, g)
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
                        me[i, j] += (N[i] * N[j]) * (g(t, x[1]) * JxW)
                    end
                end
            end
            assemble!(am, me)
        end
        return am
    end

    elit = FEIterator(Fh)
    qpit = QPIterator(Fh, (kind = :Gauss, order = 2))
    geom = geometry(Fh.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Fh), ndofs(Fh))

    integrate!(am, geom, elit, qpit, g)

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
    K = assembleK(Psih, c)
    Psi0 = fill(0.0, ndofs(Psih))
    V0 = fill(0.0, ndofs(Psih))
    for i in 1:length(locs)
        n = dofnum(Psih, 0, i, 1)
        Psi0[n] = Fic(locs[i][1])
        V0[n] = Vic(locs[i][1])
    end
    for i in bdvl
        n = dofnum(Psih, 0, i, 1)
        V0[n] = 0.0
    end
    V1 = fill(0.0, ndofs(Psih))
    L1 = fill(0.0, ndofs(Psih))
    L0 = fill(0.0, ndofs(Psih))
    # Plots
    layout = Layout(;width=700, height=700, xaxis = attr(title="x", range=[0.0, L]), yaxis = attr(title = "F", range=[-1.0, 1.0]))
    sigdig = n -> round(n*1000)/1000
    function updateplot(pl, t, xs, Psis)
        curv = scatter(;x=xs, y=Psis, mode="lines")
        plots = cat(curv; dims = 1)
        pl.plot.layout["title"] = "t = $(sigdig(t))"
        react!(pl, plots, pl.plot.layout)
        sleep(0.12)
    end
    function solve!(V1, M, K, D0, D1, L0, L1, Psi0, V0, dt, nu)
        # This assumes zero EBC
        R = (L1 + L0)/2 + (1/dt)*(M*V0) - K*(Psi0 + (dt/4)*V0) + 
            D1*((1/2)*Psi0 + (dt/4)*V0) + D0*((1/2)*Psi0)
        V1[1:nu] = ((1/dt)*M + (dt/4)*(K - D1))[1:nu, 1:nu] \ R[1:nu]
        return V1
    end
    pl = 0
    step = 0
    t = 0.0
    D0 = assembleD(Psih, g, t)
    while t < tend
        D1 = assembleD(Psih, g, t + dt)
        if mod(step, 10) == 0
            Psis = [Psi0[dofnum(Psih, 0, i, 1)] for i in 1:nshapes(ir.right)]
            if step == 0
                curv = scatter(;x=xs, y=Psis, mode="lines", name = "Sol", line_color = "rgb(155, 15, 15)")
                plots = cat(curv; dims = 1)
                pl = plot(plots, layout)
                display(pl)
            end
            updateplot(pl, t, xs, Psis)
        end
        solve!(V1, M, K, D0, D1, L0, L1, Psi0, V0, dt, nu)
        Psi1 = Psi0 + dt/2*(V1 + V0)
        t = t + dt
        step = step + 1
        copyto!(V0, V1)
        copyto!(Psi0, Psi1)
        copyto!(D0, D1)
    end
    scattersysvec!(Psih, Psi0)
    # makeattribute(Fh, "F", 1)
    # vtkwrite("gwl2-F", baseincrel(mesh), [(name = "F",)])
end

end

using .klein_gordon_var
klein_gordon_var.run()
