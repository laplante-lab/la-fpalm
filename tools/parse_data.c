/**
 * @name: John J. Ravi
 *
 *
 * TODO: implement error checking.
 *
 */

void parseImage(char* file_path, uint16 **pixels) {

  // tif files included an unknown vendor specific tag
  TIFFSetWarningHandler(NULL); // disable warnings
  TIFF* tif = TIFFOpen(file_path, "r");

  if (!tif) {
    fprintf( stderr, "%s can not be opened.\n", file_path );
    exit(1); 
  }

  uint16 h, w;
  TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &h);
  TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &w);

  if ( h != w && h != IMAGE_SIZE ) {
    fprintf( stderr, "image was: %dx%d, expected: %dx%d\n", h, w, IMAGE_SIZE, IMAGE_SIZE );
    exit(1);
  }

  for ( uint16 i = 0; i < IMAGE_SIZE; i++ ) { //height
    pixels[i] = (uint16*) malloc ( IMAGE_SIZE * sizeof( uint16 ) );
    TIFFReadScanline(tif, (void *)pixels[i], i, 0);
  }

  TIFFClose(tif);
}


