"""
    q4

Compute the solution of the Poisson equation of heat conduction with a nonzero
heat source. Quadrilateral four-node elements are used.
"""
module gwl2

using LinearAlgebra
using MeshCore.Exports
using MeshSteward.Exports
using Elfel.Exports
using PlotlyJS

L = 2.0
c = 1.0 # coefficient
N = 150;# number of subdivisions along the length of the domain
Fbc(x) = 0.0

# G(t, x) = 800*sin(5*pi*t)*(x/L)^3*(L-x)/L;# 

# Harmonic vibration
G(t, x) = 0.0; 
Fic(x) = 0.0
Vic(x) = sin(pi*x/L) + sin(3*pi*x/L)
tend = 8.0
dt = 0.01

# Uniform velocity
G(t, x) = 0.0;
Fic(x) = 0.0
Vic(x) = 1.0
tend = 8.0
dt = 0.002

# Uniform velocity, nonzero constant G
G(t, x) = -100;
Fic(x) = 0.0
Vic(x) = 1.0
tend = 8.0
dt = 0.002

# # Uniform velocity, nonzero constant G
# G(t, x) = -1000.0*(x/L)^3;
# Fic(x) = 0.0
# Vic(x) = 1.0
# tend = 8.0
# dt = 0.002

# Uniform velocity, nonzero G
# G(t, x) = +0.7*sin(0.2*t) # 
# Fic(x) = 0.0
# Vic(x) = sin(pi*x/L)
# tend = 80.0
# dt = 0.02

# Uniform velocity, nonzero constant G
G(t, x) = exp(sin(t - x));
Fic(x) = 0.0
Vic(x) = 1.0
tend = 8.0
dt = 0.002

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

function assembleD(Fh, G, t)
    function integrate!(am, geom, elit, qpit, G)
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
                        me[i, j] += (N[i] * N[j]) * (G(t, x[1]) * JxW)
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

    integrate!(am, geom, elit, qpit, G)

    return finish!(am)
end

function run()
    mesh = genmesh()
    Fh = FESpace(Float64, mesh, FEH1_L2())
    ir = baseincrel(mesh)
    bdvl = connectedv(boundary(mesh));
    locs = geometry(mesh)
    for i in bdvl
        setebc!(Fh, 0, i, 1, Fbc(locs[i]...))
    end
    numberdofs!(Fh)
    @show nu = nunknowns(Fh)
    xs = [locs[i][1] for i in 1:length(locs)]
    M = assembleM(Fh)
    K = assembleK(Fh, c)
    F0 = fill(0.0, ndofs(Fh))
    V0 = fill(0.0, ndofs(Fh))
    for i in 1:length(locs)
        n = dofnum(Fh, 0, i, 1)
        F0[n] = Fic(locs[i][1])
        V0[n] = Vic(locs[i][1])
    end
    for i in bdvl
        n = dofnum(Fh, 0, i, 1)
        V0[n] = 0.0
    end
    V1 = fill(0.0, ndofs(Fh))
    L1 = fill(0.0, ndofs(Fh))
    L0 = fill(0.0, ndofs(Fh))
    # Plots
    layout = Layout(;width=700, height=700, xaxis = attr(title="x", range=[0.0, L]), yaxis = attr(title = "F", range=[-1.0, 1.0]))
    sigdig = n -> round(n*1000)/1000
    function updateplot(pl, t, xs, Fs)
        curv = scatter(;x=xs, y=Fs, mode="lines", name = "Sol", line_color = "rgb(155, 15, 15)")
        plots = cat(curv; dims = 1)
        pl.plot.layout["title"] = "t = $(sigdig(t))"
        react!(pl, plots, pl.plot.layout)
        sleep(0.12)
    end
    function solve!(V1, M, K, D0, D1, L0, L1, F0, V0, dt, nu)
        # This assumes zero EBC
        R = (L1 + L0)/2 + (1/dt)*(M*V0) - K*(F0 + (dt/4)*V0) + 
            D1*((1/2)*F0 + (dt/4)*V0) + D0*((1/2)*F0)
        V1[1:nu] = ((1/dt)*M + (dt/4)*(K - D1))[1:nu, 1:nu] \ R[1:nu]
        return V1
    end
    pl = 0
    step = 0
    t = 0.0
    D0 = assembleD(Fh, G, t)
    while t < tend
        D1 = assembleD(Fh, G, t + dt)
        if mod(step, 10) == 0
            scattersysvec!(Fh, F0)
            makeattribute(Fh, "F", 1)
            Fa = attribute(ir.right, "F")
            Fs = [Fa[i][1] for i in 1:length(Fa)]
            if step == 0
                curv = scatter(;x=xs, y=Fs, mode="lines", name = "Sol", line_color = "rgb(155, 15, 15)")
                plots = cat(curv; dims = 1)
                pl = plot(plots, layout)
                display(pl)
            end
            updateplot(pl, t, xs, Fs)
        end
        solve!(V1, M, K, D0, D1, L0, L1, F0, V0, dt, nu)
        F1 = F0 + dt/2*(V1 + V0)
        t = t + dt
        step = step + 1
        F0, V0 = F1, V1
        D0 = D1
    end
    scattersysvec!(Fh, F0)
    # makeattribute(Fh, "F", 1)
    # vtkwrite("gwl2-F", baseincrel(mesh), [(name = "F",)])
end

end

using .gwl2
gwl2.run()
