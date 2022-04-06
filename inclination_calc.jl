function parse_value(lower, upper, name)
    println("Please enter your desired $name.")
    value = chomp(readline())
    while(true)
        try
            value = parse(typeof(lower), value)
            if (value > upper) || (value < lower)
                println("That number was out of bounds :/  The value has to be between $lower and $upper. Try again!")
                value = chomp(readline())
            else
                break
            end
        catch
            println("That was not a number you entered :/ Try again!")
            value = chomp(readline())
        end
    end
    println("Your chosen $name is $value")
    return value
end

function calc_inc_given_latitude()
    run(`clear`) 
    println("You chose calculating with your own latitude. Great choice, if I may say so!")
    lat = parse_value(-90.0, 90.0, "latitude")
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end

function choose_launchsite_and_get_latitude(launchsites)
    names = launchsites.names
    lats = launchsites.latitudes
    run(`clear`) 
    println("You chose calculating with a given launchsite. Great choice, if I may say so!")
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
    println(".....calculating inclination......")
    inc = acosd(cosd(lat)*sind(az))
    println("phew, that was hard. Your inclination is going to be $inc")
    return inc
end
function calc_az_given_launchsite(launchsites)
    while(true)
        lat = choose_launchsite_and_get_latitude(launchsites)
        inc = parse_value(-90.0,90.0,"inclination")
        if(inc<lat)
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
        else
            break
        end
    end
    az = calc_az(lat, inc)
    return az
end
function calc_az_given_latitude()
    lat = parse_value(-90.0, 90.0, "latitude")
    inc = parse_value(-90.0,90.0,"inclination")
    while(true)
        if(inc<lat)
            println("The inclination must be greater than or equal to the launch latitude, otherwise there is no mathematical solution to this problem. Please try entering different values.")
        else
            break
        end
    end
    az = calc_az(lat, inc)
    return az
end
function calc_az(lat, inc)
    println(".....calculating azimuth......")
    az = asind(cosd(inc)/cosd(lat))
    if inc!=lat
        if (cosd(inc)/cosd(lat))>0
            az2 =180- asind(cosd(inc)/cosd(lat))
            println("The problem has two possible solutions, the first Azimuth is $az, the second one is $az2")
        else
            az2 = -180-asind(cosd(inc)/cosd(lat))
            println("The problem has two possible solutions, the first Azimuth is $az, the second one is $az2")
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
    elseif option == 2
        own_lat = false
    else
        throw(BoundsError(choice,"whoa, that argument was out of bounds. how did you do that?"))
    end
    return own_lat
end
function choose_mode(launchsites)
    println("Do you want to calculate inclination or Azimuth? Type 1 for inclination, 2 for Azimuth then press enter")
    choice = parse_value(1, 2, "mode")
    #calculate inclination for given azimuth 
    if choice == 1
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_inc_given_latitude()
        else
            calc_inc_given_launchsite(launchsites)
        end
    #calculate azimuth for given inclination and latitude
    elseif choice == 2
        own_lat = latitude_or_launchsite()
        if own_lat
            calc_az_given_latitude()
        else
            calc_az_given_launchsite(launchsites)
        end
    else
        throw(BoundsError(choice,"whoa, that argument was out of bounds. how did you do that?"))
    end
end
struct Launch_sites
    names::Array
    latitudes::Array
    #cape_canaveral #28.5
    #cosmodrome  #45.9
end
function start()
    launchsites = Launch_sites([" Cape Canaveral", "Cosmodrome"],[28.5, 45.9])
    run(`clear`) 
    println("Good day, sir! Are you mayhaps in need of some inclination calculation? We have just what you need!")
    println("Here you can either calculate your inclination, given a known launch site or latitude and an azimuth or you can calculate what Azimuth you need to archieve a desired inclination.")
    choose_mode(launchsites)
end
start()
