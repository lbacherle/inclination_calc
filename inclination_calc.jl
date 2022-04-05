function parse_value(lower, upper, name)
    println("Please enter your desired $name.")
    value = chomp(readline())
    while(true)
        try
            value = parse(typeof(lower), value)
            if (value > upper) || (value < lower)
                println("That number was out of bounds :/  It should be between $lower and $upper (duh). Try again!")
                value = chomp(readline())
            else
                break
            end
        catch
            println("That was not a number you entered :/ Try again!")
            value = chomp(readline())
        end
    end
    println("Cool, your chosen $name is $value")
    return value
end

function given_latitude()
    run(`clear`) 
    println("You chose calculating with your own latitude. Great choice, if I may say so!")
    lat = parse_value(-90.0, 90.0, "latitude")
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end

function given_launchsite(launchsites)
    names = launchsites.names
    lats = launchsites.latitudes
    run(`clear`) 
    println("You chose calculating with a given launchsite. Great choice, if I may say so!")
    println("You have the following choices:\n $names \n Please type 1 for the first option, 2 for the second...you get it!")
    launchsite = parse_value(1, length(lats), "launchsite")
    lat = lats[launchsite]
    println(lats)
    println(lat)
    az = parse_value(0.0, 360.0, "azimuth")
    calc_inc(lat, az)
end
function calc_inc(lat, az)
    println(".....calculating inclination......")
    inc = acosd(cosd(lat)*sind(az))
    println("phew, that was hard. Your inclination is going to be $inc")
    return inc
end

function choose_mode(launchsites)
    println("Do you want to provide your own latitude or choose from a given launch site? Type 1 for latitude, 2 for launch site then press enter.")
    choice = parse_value(1, 2, "mode")
    if choice == 1
        given_latitude()
    elseif option == 2
        given_launchsite(launchsites)
    else
        println("whoa, you managed to parse an incompatible number. how did you do that?")
    end
end
struct Launch_sites
    names::Array
    latitudes::Array
    #cape_canaveral #28.5
    #cosmodrome  #45.9
end
function start()
    launchsites = Launch_sites([" Cape Canveral", "Cosmodrome"],[28.5, 45.9])
    run(`clear`) 
    println("Good day, sir! Are you mayhaps in need of some inclination calculation? We have just what you need!")
    println("Here you can either calculate your inclination, given a known launch site and an azimuth or for any other latitude of choice!")
    choose_mode(launchsites)
end
start()
