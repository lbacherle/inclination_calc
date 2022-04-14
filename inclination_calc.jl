"""
    inclination_calc.jl

Simple Julia program to calculate launchsite parameters.

Modes:
    1) Calculate inertial Azimuth (disregarding Earths rotation) for a launchsite or latitude and a desired inclination
    2) Calculate rotational Azimuth (taking into account Earth rotation) for a launchsite or latitude and a desired inclination
    3) Calculate resulting Inclination for a launchsite or latitude and an Azimuth angle

Conventions:
    1) All angles must be given in degrees
    2) Inclination and Latitude are defined from -90° to +90°, Azimuth from 0° to 360°
    3) Velocities must be given in m/s

Sources:
    1) Formulas: The derivation of formulas for inertial and rotational calculations can be found under https://www.orbiterwiki.org/index.php?title=Launch_Azimuth&oldid=17141
                 or https://ntrs.nasa.gov/api/citations/19980227091/downloads/19980227091.pdf ('Technical Note D-233, Determination of Azimuth Angle at Burnout for placing Satellite over a selected Earth position', T.H. Skopinski and K.G Johnson, NASA, 1960)
    2) Launchsite Latitudes: 'Space Mission Engineering - The New SMAD' by Wertz, Everett and Puschell, published 2011 by Microcosm Press

"""

using Plots
"specifies the number of calculation modes that this program currently offers"
nr_of_modes = 4

"Set desired precision of floating point values in output"
precision = 2

"Radius of Earth at equator [m]"
R_e = 6371000

"Sidereal rotational period [s]"
T_s = 86164

"struct used for storing names and latitudes of launch sites"
struct Launch_sites
    names::Array
    latitudes::Array
end

"""
    parse_value(lower, upper, name)

Asks user for an input value x, with lower <= x <= upper.

If x is not parsable to type of 'lower' or x is out of bounds the function repeats asking for input.
"""
function parse_value(lower, upper, name)
    println("Please enter your $name.")
    value = chomp(readline())
    while (true)
        try
            value = parse(typeof(lower), value)
            if (value > upper) || (value < lower)
                println("That number was out of bounds :/  The value has to be between $lower and $upper. Try again!")
                value = chomp(readline())
            else
                break
            end
        catch
            type = typeof(lower)
            println("Expected $type, but got something different. Try again!")
            value = chomp(readline())
        end
    end
    #println("Your chosen $name is $value")
    return value
end

"""
    calc_inc_given_latitude()

Asks user for input, sets values for latitude and azimuth, calls calc_inc with the given values
"""
function calc_inc_given_latitude()
    lat = parse_value(-90.0, 90.0, "latitude")
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end

function choose_launchsite_and_get_latitude(launchsites)
    names = launchsites.names
    lats = launchsites.latitudes
    println("\n##########################################################################################\n")
    println("Please choose your launchsite now.")
    println("You have the following choices:\n $names \n Type 1 for the first option, 2 for the second...")
    launchsite = parse_value(1, length(lats), "launchsite")
    lat = lats[launchsite]
    return lat
end
function calc_inc_given_launchsite(launchsites)
    lat = choose_launchsite_and_get_latitude(launchsites)
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end
function calc_inc(lat, az)
    if az <= 180
        inc = acosd(cosd(lat) * sind(az))
    else 
        inc = acosd(-cosd(lat)*sind(az))
    end
    println("The resulting inclination is ",round(inc,digits=precision),"°")
    return inc
end
function calc_rotational_az_given_launchsite(launchsites)
    lat = choose_launchsite_and_get_latitude(launchsites)
    inc = parse_value(-90.0, 90.0, "inclination")
    v = parse_value(0.0, typemax(Float32), "orbit speed")
    while (true)
        if abs(inc) < lat
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
            lat = parse_value(-90.0, 90.0, "latitude")
            inc = parse_value(-90.0, 90.0, "inclination")
            v = parse_value(0.0, typemax(Float32), "orbit speed")
        else
            break
        end
    end
    az = calc_rotational_az(lat, inc, v)
    return az
end
function calc_rotational_az_given_latitude()
    lat = parse_value(-90.0, 90.0, "latitude")
    inc = parse_value(-90.0, 90.0, "inclination")
    v = parse_value(0.0, typemax(Float32), "orbit speed")
    while (true)
        if abs(inc) < lat
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
            lat = parse_value(-90.0, 90.0, "latitude")
            inc = parse_value(-90.0, 90.0, "inclination")
            v = parse_value(0.0, typemax(Float32), "orbit speed")
        else
            break
        end
    end
    az = calc_rotational_az(lat, inc, v)
    return az
end
function calc_rotational_az(lat, inc, v)
    az_inertial = calc_az(lat, inc, false)
    v_eqrot = 465.101
    az_rot = atand((v * sind(az_inertial) - v_eqrot * cosd(inc)) / (v * cosd(az_inertial)))
    println("The rotational Azimuth needed to achieve the desired inclination of ",round(inc,digits=precision),"° is ",round(az_rot,digits=precision),"°")
    diff = abs(az_inertial - az_rot)
    println("The difference between the inertial Azimuth (",round(az_inertial,digits=precision),"°) and the rotational Azimuth (",round(az_rot,digits=precision),"°) amounts to ",round(diff,digits=precision),"°")
    return az_rot
end
function calc_az_given_latitude()
    lat = parse_value(-90.0, 90.0, "latitude")
    inc = parse_value(-90.0, 90.0, "inclination")
    while (true)
        if abs(inc) < lat
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
            lat = parse_value(-90.0, 90.0, "latitude")
            inc = parse_value(-90.0, 90.0, "inclination")
        else
            break
        end
    end
    az = calc_az(lat, inc, true)
    return az
end
function calc_az_given_launchsite(launchsites)
    lat = choose_launchsite_and_get_latitude(launchsites)
    inc = parse_value(-90.0, 90.0, "inclination")
    while (true)
        if abs(inc) < lat
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
            lat = choose_launchsite_and_get_latitude(launchsites)
            inc = parse_value(-90.0, 90.0, "inclination")
        else
            break
        end
    end
    az = calc_az(lat, inc, true)
    return az
end
function twiny(sp::Plots.Subplot)
    sp[:top_margin] = max(sp[:top_margin], 30Plots.px)
    plot!(sp.plt, inset = (sp[:subplot_index], bbox(0,0,1,1)))
    twinsp = sp.plt.subplots[end]
    twinsp[:xaxis][:mirror] = true
    twinsp[:background_color_inside] = RGBA{Float64}(0,0,0,0)
    Plots.link_axes!(sp[:yaxis], twinsp[:yaxis])
    twinsp
end
twiny(plt::Plots.Plot = current()) = twiny(plt[1])
function calc_az(lat, inc, verbose)
    az = asind(cosd(inc) / cosd(lat))
    if inc != lat
        if (cosd(inc) / cosd(lat)) > 0
            az2 = 180 - asind(cosd(inc) / cosd(lat))
            if verbose
                println("The problem has two possible solutions, the first Azimuth is ", round(az,digits=precision), "° the second one is ", round(az2,digits=precision), "° ")
            end
        else
            az2 = -180 - asind(cosd(inc) / cosd(lat))
            if verbose
                println("The problem has two possible solutions, the first Azimuth is ", round(az,digits=precision), "°, the second one is ", round(az,digits=precision), "°")
            end
        end
    else
        println("There is only one possible Azimuth that results in the given inclination, which is ", round(az,digits=precision), "°")
    end
    return az
end
function latitude_or_launchsite()
    println("Do you want to provide your own latitude or just compare different launchsites? Type 1 for latitude, 2 for launch sites.")
    choice = parse_value(1, 2, "mode")
    if choice == 1
        own_lat = true
    elseif choice == 2
        own_lat = false
    else
        throw(BoundsError(choice, "whoa, that argument was out of bounds. how did you do that?"))
    end
    return own_lat
end
function compare_propellant_mass_with_latitude(launchsites)
    lats = launchsites.latitudes
    names = launchsites.names
    #own_lat = parse_value(-90.0, 90.0, "latitude")
    #v_f = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
    #m = parse_value(0.0, typemax(Float32), "rocket mass [kg]")
    #v_e = parse_value(0.0, typemax(Float32), "exhaust velocity of gases [m/s]")
    own_lat = 0.0
    v_f = 3070
    m = 140800
    v_e = 4400
    frac = ((v_f - (2pi*R_e*cosd(own_lat))/T_s)/v_e)
    m_p=[]
    #push!(m_p,(Float64(m) * (ℯ ^frac-1)))
    m_p_own_lat = (Float64(m) * (ℯ ^frac-1))
    println("The propellant mass needed for latitude ",round(own_lat,digits=precision),"° is ",round(m_p_own_lat,digits=precision),"kg")
    println("The propellant mass for the following launchsites will also be calculated and compared in a plot: \n", names)
    for lat in launchsites.latitudes
        frac = ((v_f - (2pi*R_e*cosd(lat))/T_s)/v_e)
        push!(m_p,(Float64(m) * (ℯ ^frac-1)))
    end
    #pushfirst!(names, "input lat")
    #pushfirst!(lats, 0.0)
    p = plot(scatter(lats,m_p,mode="lines",title = "Propellant mass for different launchsite latitudes",label="p_m / launchsites"))
    xlabel!("Latitude [°]")
    ylabel!("Propellant mass [kg]")
    #twiny(plt::Plots.Plot = current()) = twiny(plt[1])
    plot!(scatter!([own_lat],[m_p_own_lat], label="p_m / given latitude",box=:on, color=:red))
    png(p,"propellant_mass_with_lat")
    display(p)
    return m_p
end

function compare_propellant_mass(launchsites)
    lats = launchsites.latitudes
    names = launchsites.names
    println("The propellant mass for the following launchsites will be calculated and compared in a plot: \n", names)
    v_f = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
    m = parse_value(0.0, typemax(Float32), "rocket mass [kg]")
    v_e = parse_value(0.0, typemax(Float32), "exhaust velocity of gases [m/s]")
    m_p=[]
    for lat in launchsites.latitudes
        frac = ((v_f - (2pi*R_e*cosd(lat))/T_s)/v_e)
        push!(m_p,(Float64(m) * (ℯ ^frac-1)))
    end
    p = plot(scatter(lats,m_p,mode="lines",title = "Propellant mass for different launchsite latitudes",label="p_m / latitudes"))
    xlabel!("Latitude [°]")
    ylabel!("Propellant mass [kg]")
    png(p,"propellant_mass")
    display(p)
    return m_p
end
function choose_mode(launchsites)
    println("\n##########################################################################################\n")
    println("Please choose your selected mode out of the $nr_of_modes modes by entering a number between 1 and $nr_of_modes.")
    choice = parse_value(1, nr_of_modes, "selection")
    #calculate inclination for given azimuth 
    if choice == 1
        println("The selected mode is calculation of the needed rotational Azimuth to achieve a certain inclination.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_rotational_az_given_latitude()
        else
            calc_rotational_az_given_launchsite(launchsites)
        end
    elseif choice == 2
        println("The selected mode is calculation of the launch Azimuth to achieve a certain inclination.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_az_given_latitude()
        else
            calc_az_given_launchsite(launchsites)
        end
    elseif choice == 3
        run(`clear`)
        println("The selected mode is calculation of the resulting inclination for a given Azimuth.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_inc_given_latitude()
        else
            calc_inc_given_launchsite(launchsites)
        end
    elseif choice == 4
        println("The selected mode is calculation of the propellant mass needed for different latitudes.")
        own_lat = latitude_or_launchsite()
        if own_lat
            compare_propellant_mass_with_latitude(launchsites)
        else
            compare_propellant_mass(launchsites)
        end
    else
        throw(BoundsError(choice, "whoa, that argument was out of bounds. how did you do that?"))
    end
end

"""
    start()

gets called at start of program. initialises launch site struct with values, prints modes and conventions on screen and calls choose_mode()
"""
function start()
    #Dict([("Cape Canaveral")])
    launchsites = Launch_sites(["Cape Canaveral", "Cosmodrome","Kennedy Space Center","Kourou Launch Center/Guyana Space Centre","LP Odissey Sea Launch","Kwajalein Missile Range","Satish Dhawan Space Centre","Jiu Quan Satellite Launch Center","Vandenberg Air Force Base"], [28.5, 45.9,28.5,5.2,0.0,9.0,13.9,40.6,34.4])
    run(`clear`)
    printstyled("LAUNCHSITE CALCULATOR\n\n", bold=true, underline=true, color=:blue)
    printstyled("Existing modes are:\n", underline=true)
    printstyled("\t1) Calculate the needed rotational Azimuth for a given latitude and inclination\n")
    printstyled("\t2) Calculate the needed inertial Azimuth for a given latitude and inclination\n")
    printstyled("\t3) Calculate the resulting inclination for a given latitude and azimuth\n")
    printstyled("\t4) Compare propellant mass needed for different latitudes\n")
    printstyled("Latitude can be provided by selecting out of a number of launchsites, or parsed directly\n")
    printstyled("\nConventions:\n", bold=true, underline=true)
    printstyled("\t-Angles must be given in degrees\n")
    printstyled("\t-Inclination and Latitude are defined from -90° to +90°, Azimuth from 0° to 360°\n")
    printstyled("\t-Velocities must be given in m/s\n")
    choose_mode(launchsites)
end
start()
