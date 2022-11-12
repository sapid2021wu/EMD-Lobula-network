function  R_output = Rect(type, input, thresthold)

R_output = input-thresthold;
if type == -1
    R_output(R_output>0) = 0;
    R_output = abs(R_output);
else 
    R_output(R_output<0) = 0;
end

