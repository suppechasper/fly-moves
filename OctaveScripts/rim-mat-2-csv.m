fid = fopen ("rims.txt");


tline = fgetl(fid);
while ischar(tline)
  load(tline);

  base = strsplit(tline, ".")(1,1){1};
  csvwrite( sprintf("%s-inner.csv", base), inner_transformed );
  csvwrite( sprintf("%s-outer.csv", base), outer_transformed );
  
  tline = fgetl(fid);
end
