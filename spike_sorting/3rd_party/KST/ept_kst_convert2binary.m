function ops = ept_kst_convert2binary(ops)

fidout = fopen(ops.fbinary, 'w');
precision = 10^round(ops.precision);
fwrite(fidout, int16(precision * ops.data), 'int16');
fclose(fidout);
ops = rmfield(ops,'data');

end