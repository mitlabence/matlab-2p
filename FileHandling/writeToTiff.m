function writeToTiff(data, path, filename)
%WRITETOTIFF Save data to a tif file.
%   data: double or uint16 raw data
%   path: ends with \, e.g. C:\tiff_images\
%   filename: e.g. example_tiff_image to save example_tiff_image.tif
output_filename = [path filename '.tif'];
disp(['Saving to ' output_filename]);
tiff_file = Tiff(output_filename, 'w8');

for i = 1:size(data,3)
        raw_data=uint16(data(:,:,i));
        tiff_file.setTag('Photometric', Tiff.Photometric.MinIsBlack);
        tiff_file.setTag('ImageWidth', size(raw_data,2));
        tiff_file.setTag('ImageLength', size(raw_data,1));
        tiff_file.setTag('PlanarConfiguration',  Tiff.PlanarConfiguration.Chunky);
        tiff_file.setTag('BitsPerSample',16);
        tiff_file.setTag('SamplesPerPixel',1);
        tiff_file.setTag('Compression', Tiff.Compression.None);
        tiff_file.write(raw_data);
        tiff_file.writeDirectory;
end
tiff_file.close();

end

