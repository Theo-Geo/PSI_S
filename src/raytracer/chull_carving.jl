function gnomonic_proj(proj,θ,φ,θmid,φmid)
    cos_Δ = sin(φmid)*sin(φ) + cos(φmid)*cos(φ)*cos(θ-θmid)
    proj[1] = R*cos(φ)*sin(θ-θmid)/cos_Δ
    proj[2] = R*(cos(φmid)*sin(φ) - sin(φmid)*cos(φ)*cos(θ-θmid))/cos_Δ
end

function carve_grid(visited,gr,source,receivers,IP)

    radius = 5
    npoints = 8
    radius_φ = deg2rad(1.0 / (111.319 * cos(gr.θ[source])) * radius)
    radius_θ = deg2rad(1.0 / 110.574 * radius)
    dΩ = 2*pi/npoints
    Ω = 0.0

    θmin, θmax = (IP.lims.lat[1]), (IP.lims.lat[2])
    φmin, φmax = (IP.lims.lon[1]), (IP.lims.lon[2])
    rmin, rmax = R + IP.lims.depth[1], R + IP.lims.depth[2]

    θmid = 0.5*(θmin+θmax)
    φmid = 0.5*(φmin+φmax)
    p = Vector{Vector{Float64}}()
    for i in 1:(1+length(receivers))*(npoints+1)
        push!(p,zeros(Float64,2))
    end

    pnode = 0
    pnode += 1
    proj = p[pnode]
    gnomonic_proj(proj,gr.θ[source],gr.φ[source],θmid,φmid)
    for receiver in receivers
        pnode += 1
        proj = p[pnode]
        gnomonic_proj(proj,gr.θ[receiver],gr.φ[receiver],θmid,φmid)
    end
    for i in 1:npoints
        pnode += 1
        proj = p[pnode]
        θ_new = gr.θ[source] + radius_θ * sin(Ω)
        φ_new = gr.φ[source] + radius_φ * cos(Ω)
        gnomonic_proj(proj,θ_new,φ_new,θmid,φmid)
        for receiver in receivers
            pnode += 1
            proj = p[pnode]
            radius_φ = deg2rad(1.0 / (111.319 * cos(gr.θ[receiver])) * radius)
            θ_new = gr.θ[receiver] + radius_θ * sin(Ω)
            φ_new = gr.φ[receiver] + radius_φ * cos(Ω)
            gnomonic_proj(proj,θ_new,φ_new,θmid,φmid)
        end
        Ω += dΩ
    end

    hull = convex_hull(p)
    vhull = VPolygon(hull)

    gproj = zeros(Float64,2)
    dθ, dφ = (θmax-θmin)/(gr.nnodes[1]-1), (φmax-φmin)/(gr.nnodes[2]-1)
    for i in 1:gr.nnodes[1], j in 1:gr.nnodes[2]
        θnode, φnode = θmin + dθ * (i-1), φmin + dφ * (j-1)
        gnomonic_proj(gproj,θnode,φnode,θmid,φmid)
        inside = (element(Singleton(gproj)) ∈ vhull) 
        if !inside 
            for k in 1:gr.nnodes[3]
                nn = LinearIndex(gr, i, j, k)
                (visited[nn] = true)
            end
        end
    end


    # xs = Float64[]
    # ys = Float64[]
    # gproj = zeros(Float64,2)
    # for i in eachindex(gr.x)
    #     if !visited[i]
    #         gnomonic_proj(gproj,gr.θ[i],gr.φ[i],θmid,φmid)
    #         push!(xs,gproj[1])
    #         push!(ys,gproj[2])
    #     end
    # end
    # h = scatter(xs,ys,markersize=1,color=:black,legend=false,aspect_ratio=:equal)
    # for i in eachindex(p)[length(receivers)+2:end]
    #     scatter!([p[i][1]],[p[i][2]],markersize=3,color=:orange,legend=false)
    # end
    # for i in eachindex(p)[2:length(receivers)+1]
    #     scatter!([p[i][1]],[p[i][2]],markersize=4,marker=:star5,color=:blue,legend=false)
    # end
    # scatter!([p[1][1]],[p[1][2]],markersize=4,marker=:utriangle,color=:red,legend=false)
    # for i in eachindex(hull)
    #     if i < length(hull)
    #         plot!([hull[i][1],hull[i+1][1]],[hull[i][2],hull[i+1][2]],c=:violet)
    #     else
    #         plot!([hull[i][1],hull[1][1]],[hull[i][2],hull[1][2]],c=:violet)
    #     end
    # end 
    # display(h)
    # readline()
end