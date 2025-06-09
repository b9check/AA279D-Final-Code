function const = getConstants(bodies)
    if ischar(bodies)
        bodies = {bodies};
    end

    const = struct();

    for i = 1:length(bodies)
        body = lower(bodies{i});
        switch body
            case 'moon'
                const.moon = utils.constants.moonConstants();
            case 'earth'
                const.earth = utils.constants.earthConstants();
            case 'mars'
                const.mars = utils.constants.marsConstants();  
            otherwise
                error('Unknown body: %s', body);
        end
    end
end
