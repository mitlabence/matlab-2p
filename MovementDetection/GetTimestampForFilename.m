function timestamp = GetTimestampForFilename()
%GETDATEFORFILENAME Get the date and time up to the precision of seconds,
%and returns it as a string. Uses format yyyymmdd_hhmise (mi = minutes, se =
%seconds)
%   The function executed on 1. February 2021, at 8:03:10 would return
%   20210201_080310

c = clock;
timestamp = strcat(num2str(c(1), '%04.f'), num2str(c(2), '%02.f'), num2str(c(3), '%02.f'), '_', num2str(c(4), '%02.f'), num2str(c(5), '%02.f'), num2str(floor(c(6)), '%02.f'));
end

