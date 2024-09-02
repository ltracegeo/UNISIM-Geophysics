function [correted_cube] = treat_border(cube)

correted_cube = cube;
for j = 1:size(cube,2)
    low_marc = false;
    high_marc = false;
    for i = 1:size(cube,1)
        
        if cube(i,j) ~= 0 
            if low_marc==false;
                correted_cube(1:i,j) = cube(i+1,j);
                low_marc = true;
            end
        end
        
        if (cube(i,j) == 0) || isnan(cube(i,j))
            if high_marc==false && low_marc==true;
                correted_cube(i:end,j) = cube(i-1,j);
                high_marc = true;
            end
            
        end
    end
end
media = 10000;

for j = 1:size(correted_cube,2)

    for i = 1:size(correted_cube,1)
        if correted_cube(i,j) == 0 
            correted_cube(i,j) = media;
        end
    end
    
end