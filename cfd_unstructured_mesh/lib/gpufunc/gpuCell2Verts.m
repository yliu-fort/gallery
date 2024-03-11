function output = gpuCell2Verts( field,nodeInterpOperator)
% high performance interpolation function
% interp cell value to nodes


    output = nodeInterpOperator*field;
    
end
