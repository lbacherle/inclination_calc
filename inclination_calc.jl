###############
nr_of_modes = 3
#set desired precision of floating point values
precision = 2
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
    inc = acosd(cosd(lat) * sind(az))
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
    println(".....calculating azimuth......")
    az_inertial = calc_az(lat, inc, false)
    v_eqrot = 465.101
    az_rot = atand((v * sind(az_inertial) - v_eqrot * cosd(inc)) / (v * cosd(az_inertial)))
    println("The rotational Azimuth needed to achieve the desired inclination of ",round(inc,digits=precision),"° is ",round(az_rot,digits=precision),"°")
    diff = abs(az_inertial - az_rot)
    println("The difference between the inertial Azimuth (",round(az_inertial,digits=precision),"°) and the rotational Azimuth (",round(az_rot,digits=precision),"°) is ",round(diff,digits=precision),"°")
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
function calc_az(lat, inc, verbose)
    if verbose
        println(".....calculating azimuth......")
        sleep(1)
    end
    az = asind(cosd(inc) / cosd(lat))
    if inc != lat
        if (cosd(inc) / cosd(lat)) > 0
            az2 = 180 - asind(cosd(inc) / cosd(lat))
            if verbose
                println("The problem has two possible solutions, the first Azimuth is $az, the second one is $az2")
            end
        else
            az2 = -180 - asind(cosd(inc) / cosd(lat))
            if verbose
                println("The problem has two possible solutions, the first Azimuth is $az, the second one is $az2")
            end
        end
    else
        println("There is only one possible Azimuth that results in the given inclination, which is $az")
    end
    return az
end
function latitude_or_launchsite()
    println("Do you want to provide your own latitude or choose from a given launch site? Type 1 for latitude, 2 for launch site.")
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
function choose_mode(launchsites)
    println("\n##########################################################################################\n")
    println("Please choose your selected mode out of the $nr_of_modes modes.")
    choice = parse_value(1, 3, "mode")
    #calculate inclination for given azimuth 
    if choice == 3
        run(`clear`)
        println("The selected mode is calculation of the resulting inclination for a given Azimuth.")
        sleep(1)
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_inc_given_latitude()
        else
            calc_inc_given_launchsite(launchsites)
        end
    elseif choice == 2
        println("The selected mode is calculation of the launch Azimuth to achieve a certain inclination.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_az_given_latitude()
        else
            calc_az_given_launchsite(launchsites)
        end
    elseif choice == 1
        println("The selected mode is calculation of the needed rotational Azimuth to achieve a certain inclination.")
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_rotational_az_given_latitude()
        else
            calc_rotational_az_given_launchsite(launchsites)
        end
    else
        throw(BoundsError(choice, "whoa, that argument was out of bounds. how did you do that?"))
    end
end
struct Launch_sites
    names::Array
    latitudes::Array
    #cape_canaveral #28.5
    #cosmodrome  #45.9
end
function start()
    launchsites = Launch_sites(["Cape Canaveral", "Cosmodrome"], [28.5, 45.9])
    run(`clear`)
    printstyled("LAUNCHSITE CALCULATOR\n\n", bold=true, underline=true, blink=true, color=:blue)
    printstyled("Existing modes are:\n", underline=true)
    printstyled("\t1) Calculate the needed rotational Azimuth for a given latitude and inclination\n")
    printstyled("\t2) Calculate the needed inertial Azimuth for a given latitude and inclination\n")
    printstyled("\t3) Calculate the resulting inclination for a given latitude and azimuth\n")
    printstyled("Latitude can be provided by selecting out of a number of launchsites, or parsed directly\n")
    printstyled("\nConventions:\n", bold=true, underline=true)
    printstyled("\t-Angles must be given in degrees\n")
    printstyled("\t-Inclination and Latitude is defined from -90° to +90°, Azimuth from 0° to 360°\n")
    choose_mode(launchsites)
end
start()
