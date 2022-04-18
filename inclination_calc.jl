"""
    inclination_calc.jl

Simple Julia program to calculate launchsite parameters.

Modes:
    1) Calculate rotational Azimuth (taking into account Earth rotation) for a launchsite or latitude and a desired Inclination
    2) Calculate inertial Azimuth (disregarding Earths rotation) for a launchsite or latitude and a desired inclination
    3) Calculate resulting Inclination for a launchsite or latitude and an Azimuth angle
    4) Compare neccessary propellant mass as a function of latitude (assuming an eastwards launch)

Input parameters:
    Mode 1: Latitude/Launchsite selection (-90.0° -> +90°), Inclination (-90.0° -> +90.0°), Final Orbit Velocity (0 m/s -> max(Float32))
    Mode 2: Latitude/Launchsite selection (-90.0° -> +90°), Inclination (-90.0° -> +90.0°)
    Mode 3: Latitude/Launchsite selection (-90.0° -> +90°), Azimuth (0° -> 360°)
    Mode 4: Latitude/Launchsite selection (-90.0° -> +90°), Final Orbit Velocity (0 m/s -> max(Float32)), Rocket mass (0kg -> max(Float32)), Exhaust Velocity (0 m/s -> max(Float32))

Input:
    If file called 'values.csv' in the correct format (see example values.csv) exists, reads the values from this csv and prints all output to files.
    If this file does not exist, starts an interactive mode and allows user to freely input variables.

Conventions:
    1) All angles must be given in degrees
    2) Inclination and Latitude are defined from -90° to +90°, Azimuth from 0° to 360°
    3) Velocities must be given in m/s
    
Dependencies:
    1) The Plots package is needed to show the resulting plots. Installation can be done by typing julia into the console to start the REPL, then typing the following:
        julia> using package
        julia> Pkg.add("Plots")

Sources:
    1) Formulas: The derivation of formulas for inertial and rotational calculations can be found under https://www.orbiterwiki.org/index.php?title=Launch_Azimuth&oldid=17141
                 or https://ntrs.nasa.gov/api/citations/19980227091/downloads/19980227091.pdf ('Technical Note D-233, Determination of Azimuth Angle at Burnout for placing Satellite over a selected Earth position', T.H. Skopinski and K.G Johnson, NASA, 1960)
    2) Launchsite Latitudes: 'Space Mission Engineering - The New SMAD' by Wertz, Everett and Puschell, published 2011 by Microcosm Press

Additional Information:
    Developer:       Lisa Bacherle
    Version:         1.0
    Date of release: 17.04.2022
    Developed using: Julia 1.7.2 (2022-02-06)
"""

using Plots
using CSV
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

If x is not parsable to type of 'lower' or x is out of bounds the function repeats asking for input until both conditions are satisfied.
Returns x.
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
    choose_launchsite_and_get_latitude(launchsites)
    
Prints launchsite names, asks user for input on which launchsite he wants to choose. 
Returns the latitude of the chosen launchsite.
"""
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

"""
    latitude_or_launchsite()

Asks if user wants to provide own latitude or choose from a list of given launchsites. Returns true for first option, false for second.
"""
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

"""
    calc_inc(lat, az)

Calculates the resulting inclination of the two arguments. Prints and returns the inclination.
"""
function calc_inc(lat, az)
    if az <= 180
        inc = acosd(cosd(lat) * sind(az))
    else 
        inc = acosd(-cosd(lat)*sind(az))
    end
    println("The resulting inclination is ",round(inc,digits=precision),"°")
    return inc
end

"""
    calc_inc(lat, az,out)

Calculates the resulting inclination of the two arguments. Write to out file and returns the inclination.
"""
function calc_inc(lat, az,out)
    if az <= 180
        inc = acosd(cosd(lat) * sind(az))
    else 
        inc = acosd(-cosd(lat)*sind(az))
    end
    write(out,"The resulting inclination for the latitude ",repr(lat)," degrees and the Azimuth ", repr(az)," degrees is ",repr(round(inc,digits=precision))," degrees")
    return inc
end

"""
    calc_inc_given_latitude()

Asks user for input, sets values for latitude and azimuth, calls calc_inc with the given values.
"""
function calc_inc_given_latitude()
    lat = parse_value(-90.0, 90.0, "latitude")
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end

"""
    calc_inc_given_launchsite

Asks for user input for launchsite and azimuth. Calls calc_inc with these values.
"""
function calc_inc_given_launchsite(launchsites)
    lat = choose_launchsite_and_get_latitude(launchsites)
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end

"""
    calc_rotational_az(lat, inc, v)

Assumes that calculation of a rotational Azimuth is possible for the passed values. Calculates inertial azimuth using calc_az and uses the result to calculate rotational azimuth.
Prints and returns rotational Azimuth.
"""
function calc_rotational_az(lat, inc, v)
    az_inertial = calc_az(lat, inc, false)
    v_eqrot = 465101.0
    az_rot = atand((v * sind(az_inertial) - v_eqrot * cosd(inc)) / (v * cosd(az_inertial)))
    println("The rotational Azimuth needed to achieve the desired inclination of ",round(inc,digits=precision),"° equals ",round(az_rot,digits=precision),"°")
    diff = abs(az_inertial - az_rot)
    println("The difference between the inertial Azimuth (",round(az_inertial,digits=precision),"°) and the rotational Azimuth (",round(az_rot,digits=precision),"°) equals ",round(diff,digits=precision),"°")
    return az_rot
end

"""
    calc_rotational_az(lat, inc, v, out)

Assumes that calculation of a rotational Azimuth is possible for the passed values. Calculates inertial azimuth using calc_az and uses the result to calculate rotational azimuth.
Writes to file and returns rotational Azimuth.
"""
function calc_rotational_az(lat, inc, v, out)
    az_inertial = calc_az(lat, inc, false)
    v_eqrot = 465101.0
    az_rot = atand((v * sind(az_inertial) - v_eqrot * cosd(inc)) / (v * cosd(az_inertial)))
    write(out,"\nThe rotational Azimuth needed to achieve the desired inclination of ")
    show(out,round(inc,digits=precision))
    write(out,"° equals ")
    show(out,repr(round(az_rot,digits=precision)))
    write(out,"degrees\n")
    diff = abs(az_inertial - az_rot)
    write(out,"The difference between the inertial Azimuth (",repr(round(az_inertial,digits=precision))," degrees) and the rotational Azimuth (",repr(round(az_rot,digits=precision))," degrees) equals ",repr(round(diff,digits=precision))," degrees\n")
    return az_rot
end

"""
    calc_rotational_az_given_launchsite(launchsites)

Asks for user input for launchsite, inclination and velocity. Checks that inclination is bigger or equal than latitude, if not keeps asking for input until condition is fulfilled.
Calls calc_rotational_az with the parsed values.
"""
function calc_rotational_az_given_launchsite(launchsites)
    lat = choose_launchsite_and_get_latitude(launchsites)
    inc = parse_value(-90.0, 90.0, "inclination")
    v = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
    while (true)
        if abs(inc) < lat
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
            lat = choose_launchsite_and_get_latitude(launchsites)
            inc = parse_value(-90.0, 90.0, "inclination")
            v = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
        else
            break
        end
    end
    az = calc_rotational_az(lat, inc, v)
    return az
end

"""
    calc_rotational_az_given_latitude()

Asks for user input for latitude, inclination and velocity. Checks that inclination is bigger or equal than latitude, if not keeps asking for input until condition is fulfilled.
Calls calc_rotational_az with the parsed values.
"""
function calc_rotational_az_given_latitude()
    lat = parse_value(-90.0, 90.0, "latitude")
    inc = parse_value(-90.0, 90.0, "inclination")
    v = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
    while (true)
        if abs(inc) < lat
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
            lat = parse_value(-90.0, 90.0, "latitude")
            inc = parse_value(-90.0, 90.0, "inclination")
            v = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
        else
            break
        end
    end
    az = calc_rotational_az(lat, inc, v)
    return az
end

"""
    calc_az(lat,inc,verbose)

Assumes that calculation of inertial Azimuth is possible for the passed parameters. Calculates one or, if possible, two solutions for inertial Azimuth.
If verbose, prints them and returns the first solution.
If verbose is not true, returns the first solution to the Azimuth, which corresponds to an launch headed east.
"""
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

"""
    calc_az_given_latitude()

Asks user for input on and sets values for latitude and inclination, calls calc_az with the given values.
"""
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

"""
    calc_az_given_launchsite(launchsites)

Asks user for input on and sets values for launchsite and inclination, calls calc_az with the given values.
"""
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

"""
    compare_propellant_mass_with_latitude(launchsites)

Asks for user input on latitude, final orbit velocity, rocket mass and exhaust velocity. Calculates the propellant mass neccessary to achieve the desired orbit velocity.
Does the same for every launchsite in launchsites. Plots the results and saves the plot to a file to facilitate comparison.
Returns the necessary propellant mass for the desired orbit velocity.
"""
function compare_propellant_mass_with_latitude(launchsites)
    lats = launchsites.latitudes
    names = launchsites.names

    #input values
    own_lat = parse_value(-90.0, 90.0, "latitude")
    v_f = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
    m = parse_value(0.0, typemax(Float32), "rocket mass [kg]")
    v_e = parse_value(0.0, typemax(Float32), "exhaust velocity of gases [m/s]")

    #calculation of p_m for input latitude
    frac = ((v_f - (2pi*R_e*cosd(own_lat))/T_s)/v_e)
    m_p=[]
    push!(m_p,(Float64(m) * (ℯ ^frac-1)))
    pushfirst!(names, "input lat")
    pushfirst!(lats, 0.0)
    println("The propellant mass needed for latitude ",round(own_lat,digits=precision),"° equals ",round(m_p[1],digits=precision),"kg")
    println("The propellant mass for the following launchsites will also be calculated and compared in a plot: \n", names)

    #calculation of p_m for launchsites
    for lat in launchsites.latitudes
        frac = ((v_f - (2pi*R_e*cosd(lat))/T_s)/v_e)
        push!(m_p,(Float64(m) * (ℯ ^frac-1)))
    end

    #plotting
    p = plot(scatter(lats,m_p,mode="lines",title = "propellant mass / latitudes (eastwards launch)",label="p_m / launchsites"))
    xlabel!("Latitude [°]")
    ylabel!("Propellant mass p_m [kg]")
    plot!(scatter!([own_lat],[m_p[1]], label="p_m / chosen latitude",box=:on, color=:red))
    png(p,"propellant_mass_with_lat")
    display(p)

    return m_p[1]
end

"""
    compare_propellant_mass(launchsites)

Asks for user input on final orbit velocity, rocket mass and exhaust velocity. For every launchsite in launchsites, calculates the propellant mass neccessary to achieve the desired orbit velocity.
Plots the results and saves the plot to a file to facilitate comparison.
Returns an array containing the propellant mass for every launchsite.
"""
function compare_propellant_mass(launchsites)
    lats = launchsites.latitudes
    names = launchsites.names
    m_p=[]
    println("The propellant mass for the following launchsites will be calculated and compared in a plot: \n", names)

    #input values
    v_f = parse_value(0.0, typemax(Float32), "final orbit velocity [m/s]")
    m = parse_value(0.0, typemax(Float32), "rocket mass [kg]")
    v_e = parse_value(0.0, typemax(Float32), "exhaust velocity of gases [m/s]")
    
    #calculation of m_p for launchsites
    for lat in launchsites.latitudes
        frac = ((v_f - (2pi*R_e*cosd(lat))/T_s)/v_e)
        push!(m_p,(Float64(m) * (ℯ ^frac-1)))
    end

    #plotting
    p = plot(scatter(lats,m_p,mode="lines",title = "propellant mass / latitudes (eastwards launch)",label="p_m / latitudes"))
    xlabel!("Latitude [°]")
    ylabel!("Propellant mass p_m [kg]")
    png(p,"propellant_mass")
    display(p)

    return m_p
end

"""
    compare_propellant_mass(launchsites, values)

Noninteractive version of compare_propellant_mass.
Gets from values the following: final orbit velocity, rocket mass and exhaust velocity. For every launchsite in launchsites, calculates the propellant mass neccessary to achieve the desired orbit velocity.
Also calculates the difference in necessary propellant mass for every launchsite compared to the minimum propellant mass needed at equator.
Plots the results and saves the plot to a file to facilitate comparison.
Returns an array containing the propellant mass for every launchsite.
"""
function compare_propellant_mass(launchsites, values)
    lats = launchsites.latitudes
    names = launchsites.names
    m_p=[]
    delta_mp = []
    println("The propellant mass for the following launchsites will be calculated and compared in a plot: \n", names)
    pushfirst!(lats, 0)
    #input values
    v_f = values[4]
    m = values[5]
    v_e = values[6]
    frac = ((v_f - (2pi*R_e*cosd(0))/T_s)/v_e)
    mp_0 = Float64(m) * (ℯ ^frac-1)
    #calculation of m_p for launchsites
    for lat in launchsites.latitudes
        frac = ((v_f - (2pi*R_e*cosd(lat))/T_s)/v_e)
        push!(m_p,(Float64(m) * (ℯ ^frac-1)))
        #calculate difference between propellant mass needed at launchsite and at equator
        push!(delta_mp,abs((Float64(m) * (ℯ ^frac-1)-mp_0)))
    end

    #plotting
    p = plot(scatter(lats,m_p,mode="lines",title = "propellant mass / latitudes (eastwards launch)",label="p_m / latitudes"))
    xlabel!("Latitude [°]")
    ylabel!("Propellant mass p_m [kg]")
    png(p,"propellant_mass_launchsites")
    display(p)

    #plotting plot of delta m_p
    p2 = plot(scatter(lats,delta_mp,mode="lines",title = "propellant mass / latitudes (eastwards launch)",label="Δp_m / latitudes"))
    xlabel!("Latitude [°]")
    ylabel!("Difference in propellant mass Δp_m [kg]")
    png(p2,"delta_propellant_mass_launchsites")
    display(p2)
    return m_p
end

"""
    compare_propellant_mass_all_latitudes(values)

Gets from values the following: final orbit velocity, rocket mass and exhaust velocity. For every latitude from 0° to 90°, calculates the propellant mass neccessary to achieve the desired orbit velocity.
Also calculates the difference in necessary propellant mass for every launchsite compared to the minimum propellant mass needed at equator.
Plots the results and saves plot to a file to facilitate comparison.
Returns an array containing the propellant mass for every latitude.
"""
function compare_propellant_mass_all_latitudes(values)
    m_p=[]
    delta_mp = []
    lats = []
    #input values
    v_f = values[4]
    m = values[5]
    v_e = values[6]
    frac = ((v_f - (2pi*R_e*cosd(0))/T_s)/v_e)
    mp_0 = Float64(m) * (ℯ ^frac-1)
    #calculation of m_p for launchsites
    for lat in 0:90
        push!(lats,lat)
        frac = ((v_f - (2pi*R_e*cosd(lat))/T_s)/v_e)
        push!(m_p,(Float64(m) * (ℯ ^frac-1)))
        #calculate difference between propellant mass needed at launchsite and at equator
        push!(delta_mp,abs((Float64(m) * (ℯ ^frac-1)-mp_0)))
    end

    #plotting
    p = plot(lats,m_p,title = "propellant mass / latitudes (eastwards launch)",label="p_m / latitudes")
    xlabel!("Latitude [°]")
    ylabel!("Propellant mass p_m [kg]")
    png(p,"propellant_mass_all_latitudes")
    display(p)

    #plotting plot of delta m_p
    p2 = plot(lats,delta_mp,title = "propellant mass / latitudes (eastwards launch)",label="Δp_m / latitudes")
    xlabel!("Latitude [°]")
    ylabel!("Difference in propellant mass Δp_m [kg]")
    png(p2,"delta_propellant_mass_all_latitudes")
    display(p2)
    return m_p
end

"""
    choose_mode(launchsites)

Asks user for input on which of the 4 modes to use. Asks whether user wants to provide own latitude or choose from the available launchsites.
Calls the appropriate function for the chosen mode.
"""
function choose_mode(launchsites)
    println("\n##########################################################################################\n")
    println("Please choose your selected mode out of the $nr_of_modes modes by entering a number between 1 and $nr_of_modes.")
    choice = parse_value(1, nr_of_modes, "selection")

    if choice == 1
        println("The selected mode is calculation of the rotational Azimuth neccessary to achieve a certain inclination.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_rotational_az_given_latitude()
        else
            calc_rotational_az_given_launchsite(launchsites)
        end
    elseif choice == 2
        println("The selected mode is calculation of the Azimuth neccessary to achieve a certain inclination.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_az_given_latitude()
        else
            calc_az_given_launchsite(launchsites)
        end
    elseif choice == 3
        println("The selected mode is calculation of the resulting inclination for a given Azimuth.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_inc_given_latitude()
        else
            calc_inc_given_launchsite(launchsites)
        end
    elseif choice == 4
        println("The selected mode is calculation of the propellant mass needed as a function of latitude")
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
    read_from_csv(file)

Reads values from a csv and converts them to a matrix. Returns the matrix.
"""
function read_from_csv(file)
    values = file |> CSV.Tables.matrix
    return values
end

function init_interactive(launchsites)
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

"""
    init_noninteractive(launchsites, file, out)

Initialises values for noninteractive mode by reading them from file.
Calls all functions to calculate the 4 modes.
Closes files.
"""
function init_noninteractive(launchsites, file, out)
    write(out,"LAUNCHSITE CALCULATOR\n\n")
    printstyled("LAUNCHSITE CALCULATOR\n\n", bold=true, underline=true, color=:blue)
    printstyled("The following things will be calculated and the results written to a file:\n", underline=true)
    printstyled("\t1) Calculating necessary rotational Azimuth for a given latitude and inclination\n")
    printstyled("\t2) Calculating necessary inertial Azimuth for a given latitude and inclination\n")
    printstyled("\t3) Calculating resulting inclination for a given latitude and azimuth\n")
    printstyled("\t4) Comparing propellant mass needed for different latitudes\n")
    all_values = read_from_csv(file)
    calc_az(all_values[2,1],all_values[2,2],true)
    calc_rotational_az(all_values[1,1],all_values[1,2],all_values[1,4],out)
    calc_inc(all_values[3,1],all_values[3,3],out)
    compare_propellant_mass(launchsites,all_values[4,:])
    compare_propellant_mass_all_latitudes(all_values[4,:])
    close(out)
end
"""
    start()

gets called at start of program. Initialises launch site struct with values.
Tries to find and open a file called values.csv. If this file exists, generates an output file and starts noninteractive mode by calling init_noninteractive.
If this file does not exist, starts interactive mode by calling init_interactive().
"""
function start()
    launchsites = Launch_sites(["Cape Canaveral", "Cosmodrome","Kennedy Space Center","Kourou Launch Center/Guyana Space Centre","LP Odissey Sea Launch","Kwajalein Missile Range","Satish Dhawan Space Centre","Jiu Quan Satellite Launch Center","Vandenberg Air Force Base","Santa Maria Island, Azores"], [28.5, 45.9,28.5,5.2,0.0,9.0,13.9,40.6,34.4,37.0])
    run(`clear`)
    try
        file = CSV.File(open("values.csv"))
        touch("output_values.txt")
        out = open("output_values.txt","w")
        init_noninteractive(launchsites, file,out)
    catch
        init_interactive(launchsites)
    end
end

start()
