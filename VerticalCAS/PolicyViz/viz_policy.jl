export viz_policy

function viz_policy(nnet_folder::AbstractString,table_folder::AbstractString="",vnum::Int64=2, hu::Int64=45, epoch::Int64=1200)
    
    batch_size=500
    pra  = 1
    currentPra = 1
    
    if epoch>0
        nnet = NNet(@sprintf("%s/VerticalCAS_r15_pra%02d_v%d_%dHU_%04d.nnet",nnet_folder,pra,vnum, hu, epoch))
    else
        nnet = NNet(@sprintf("%s/VerticalCAS_r15_pra%02d_v%d_%dHU.nnet",nnet_folder,pra,vnum, hu))
    end

    grid = RectangleGrid(relhs,dh0s,dh1s,taus)
    policy = []
    if table_folder != ""
        alpha = h5open(@sprintf("%s/VerticalCAS_r15_TrainingData_%02d.h5", table_folder,pra), "r") do file
            read(file, "y")
        end
        alpha = alpha';
        policy = read_policy(NEXT_ACTIONS[pra], alpha)
    end
    
    ra_1 = RGB(1.,1.,1.) # white
    ra_2 = RGB(.0,.0,.5) # navy
    ra_3 = RGB(.0,.600,.0) # dark green
    ra_4 = RGB(.5,.5,.5) # grey
    ra_5 = RGB(.7,.9,.0) # neon green
    ra_6 = RGB(1.0,0.56,0.0) # orange
    ra_7 = RGB(1.0,0.72,0.76) # pink
    ra_8 = RGB(.9,.500,.9) # magenta
    ra_9 = RGB(0.9,.9,.5) # yellow
    ra_10 = RGB(0.2,.7,0.9) # neon blue
    ra_11 = RGB(0.6,0.0,0.0) # red
    colors = [ra_1;ra_2;ra_3;ra_4;ra_5;ra_6;ra_7;ra_8;ra_9;ra_10;ra_11]
    
    hp = 100
    
    @manipulate for nbin = 100,
        savePlot = [false,true],
        xmin = 6.0,
        xmax = 40.0,
        ymin = -1000.0,
        ymax = 1000.0,
        dh0 = 0.0, 
        dh1 = 0.0, 
        pra = pas,
        w = -1.0,
        vlo = 0.0,
        showRegion = [false,true]
        
        
        if pra != currentPra
            if table_folder !=""
                alpha = h5open(@sprintf("%s/VerticalCAS_r15_TrainingData_%02d.h5",table_folder,pra), "r") do file
                    read(file, "y")
                end
                alpha = alpha';
                policy = read_policy(NEXT_ACTIONS[pra], alpha)
            end
            
            if epoch>0
                nnet = NNet(@sprintf("%s/VerticalCAS_r15_pra%02d_v%d_%dHU_%04d.nnet",nnet_folder,pra,vnum,hu,epoch))
            else
                nnet = NNet(@sprintf("%s/VerticalCAS_r15_pra%02d_v%d_%dHU.nnet",nnet_folder,pra,vnum, hu))
            end
            currentPra = pra
        end
        
        
        #Load table with the inputs needed to plot the heat map
        inputsNet = []
        inputsNet= zeros(nbin*nbin,num_inputs(nnet)) 
        ind = 1
        for i=linspace(xmin, xmax, nbin)
            for j=linspace(ymin, ymax, nbin)
                tau = i
                relh = j
                if num_inputs(nnet) == 4
                    inputsNet[ind,:] = [relh,dh0,dh1,tau];
                elseif num_inputs(nnet) == 5
                    inputsNet[ind,:] = [relh,dh0,dh1,ra,tau];
                end
                ind = ind+1
            end
        end
        
        #Calculate all of the Q values from the input array
        q_nnet = zeros(nbin*nbin,num_outputs(nnet));
        ind = 1
        
        while ind+batch_size<nbin*nbin
            input = inputsNet[ind:(ind+batch_size-1),:]'
            output = evaluate_network_multiple(nnet,input) 
            q_nnet = [q_nnet[1:(ind-1),:];output';q_nnet[ind+batch_size:end,:]]
            ind=ind+batch_size
        end
        input = inputsNet[ind:end,:]'
        output = evaluate_network_multiple(nnet,input)
        q_nnet = [q_nnet[1:(ind-1),:];output']

        
        # Q Table Heat Map
        function get_heat1(x::Float64, y::Float64)
            tau = x
            relh = y
            qvals = evaluate(policy, get_belief([relh,dh0,dh1,tau],grid,false))
            return NEXT_ACTIONS[pra][indmax(qvals)]
        end # function get_heat1
        
        ind = 1
        #Neural Net Heat Map
        function get_heat2(x::Float64, y::Float64)                        
            qvals = q_nnet[ind,:]
            ind +=1
            return NEXT_ACTIONS[pra][indmax(qvals)]
        end # function get_heat2
        
        
        g = GroupPlot(2, 1, groupStyle = "horizontal sep=1cm, vertical sep=2.5cm")
        if table_folder!=""
            g = GroupPlot(2, 2, groupStyle = "horizontal sep=4cm, vertical sep=2.5cm")
            if showRegion
                spacing = (ymax-ymin)/30.0
                push!(g, Axis(vcat([
                Plots.Image(get_heat1, (xmin, xmax), (ymin, ymax), zmin = 1, zmax = 11, xbins = nbin, ybins = nbin, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false)],[
                Plots.Linear([xmin,xmax],vlo*[xmin, xmax]-w*hp + w*spacing*i, style="solid, black, no marks") for i=0:60],[
                Plots.Linear([xmin,xmax],vlo*[xmin, xmax]-w*hp, style="solid, black, thick, no marks")]), 
                xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax,width="10cm", height="8cm", xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Original Table Advisories"))
            else
                push!(g, Axis([
                    Plots.Image(get_heat1, (xmin, xmax), 
                    (ymin, ymax), 
                    zmin = 1, zmax = 11,
                    xbins = nbin, ybins = nbin,
                    colormap = ColorMaps.RGBArrayMap(colors), colorbar=false),
                    ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax, width="10cm", height="8cm", xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Original Table Advisories"))
            end
        end
        if showRegion
            spacing = (ymax-ymin)/30.0
                push!(g, Axis(vcat([
                Plots.Image(get_heat2, (xmin, xmax), (ymin, ymax), zmin = 1, zmax = 11, xbins = nbin, ybins = nbin, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false)],[
                Plots.Linear([xmin,xmax],vlo*[xmin, xmax]-w*hp + w*spacing*i, style="solid, black, no marks") for i=0:60],[
                Plots.Linear([xmin,xmax],vlo*[xmin, xmax]-w*hp, style="solid, black, thick, no marks")]), 
                xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax,width="10cm", height="8cm", xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Neural Network Advisories"))
        else
            push!(g, Axis([
                Plots.Image(get_heat2, (xmin, xmax), (ymin, ymax), zmin = 1, zmax = 11, xbins = nbin, ybins = nbin, colormap = ColorMaps.RGBArrayMap(colors), colorbar=false), 
                ], xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax,width="10cm", height="8cm", xlabel="Tau (s)", ylabel="Relative Alt (ft)", title="Neural Network Advisories"))
        end
        if savePlot
            PGFPlots.save("PolicyPlot.tex", g, include_preamble=true)
        else
            bg_colors = [RGB(1.0,1.0,1.0)]
            define_color("ra_1",  [0.95,0.95,0.95])
            define_color("ra_2",  [0.,0.,0.5])
            define_color("ra_3",  [.0,.600,.0])
            define_color("ra_4",  [.5,.5,.5])
            define_color("ra_5",  [.7,.9,.0])
            define_color("ra_6",  [1.0,0.56,0.0])
            define_color("ra_7",  [1.0,0.72,0.76])
            define_color("ra_8",  [.9,.500,.9])
            define_color("ra_9",  [0.9,.9,.5])
            define_color("ra_10", [0.2,.7,0.9])
            define_color("ra_11", [0.6,0.0,0.0])
            xx = [-1.5,-1.5,-1.5, -1.5, -1.5, -1.5, -1.5, 0.4 ,.4, .4, .4]
            yy = [1.65,1.15,0.65, 0.15, -0.35, -0.85, -1.35, 1.65, 1.15, 0.65, 0.15]
            zz = ["ra_1","ra_2","ra_3","ra_4","ra_5","ra_6","ra_7","ra_8","ra_9","ra_10","ra_11"]
            sc = string("{ra_1={mark=square, style={black, mark options={fill=ra_1}, mark size=6}},",
                         "ra_2={style={ra_2, mark size=6}},",
                         "ra_3={style={ra_3, mark size=6}},",
                         "ra_4={style={ra_4, mark size=6}},",
                         "ra_5={style={ra_5, mark size=6}},",
                         "ra_6={style={ra_6, mark size=6}},",
                         "ra_7={style={ra_7, mark size=6}},",
                         "ra_8={style={ra_8, mark size=6}},",
                         "ra_9={style={ra_9, mark size=6}},",
                         "ra_10={style={ra_10, mark size=6}},",
                         "ra_11={style={ra_11, mark size=6}}}"
                )

            f = (x,y)->x*exp(-x^2-y^2)
            push!(g, Axis([
                Plots.Image(f, (-2,2), (-2,2),colormap = ColorMaps.RGBArrayMap(bg_colors),colorbar=false),
                Plots.Scatter(xx, yy, zz, scatterClasses=sc),
                Plots.Node("RA 1: COC ",0.15,0.915,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 2: DNC-1/4 ",0.15,0.790,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 3: DND-1/4 ",0.15,0.665,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 4: Maintain-1/4 ",0.15,0.540,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 5: D15000-1/4 ",0.15,0.415,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 6: C1500-1/4 ",0.15,0.290,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 7: D1500-1/3 ",0.15,0.165,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 8:  C1500-1/3 ",0.63,0.915,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 9:  D2500-1/3 ",0.63,0.790,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 10: C2500-1/3 ",0.63,0.665,style="black,anchor=west", axis="axis description cs"),
                Plots.Node("RA 11: MTLO-1/4 ",0.63,0.540,style="black,anchor=west", axis="axis description cs"),
                ],width="10cm",height="8cm", hideAxis =true, title="KEY"))
        end
    end # for p_int, v0, v1, pa, ta
end # function viz_pairwise_policy