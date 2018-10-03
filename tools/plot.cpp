#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>

#include "plot.h"
#include <pthread.h>

#include <string>

using namespace std;

void gnuprint(FILE *gp, float *x, int N) {     
  for (int i=0; i<N; i++)
    fprintf(gp, "%g\n", x[i]);

  fflush(gp);
  fprintf(gp, "e\n");
}

void gnuprint(FILE *gp, float *x, float *y, float *z, int N) { 
  for (int i=0; i<N; i++)
    fprintf(gp, "%g %g %g\n", x[i], y[i], z[i]);

  fflush(gp);
  fprintf(gp, "e\n");
}

float cal_mean(int n, float *x) {
  float sum = 0;

  for(int i=0; i<n; i++)
    sum+=x[i];

  return((float)sum/n);
}

void *plotFrequencyEmissionsBg( void *data ) {

  struct fe_data *func_data = (struct fe_data *) data;

  FILE *f = popen("gnuplot", "w");

  fprintf(f, "set terminal png small\n" 
              //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
             "set output '%s'\n"
             "unset key\n"
             "set multiplot layout 1,2\n"
             "binwidth=1\n" 
             "bin(x,width)=width*floor(x/width)\n" 
             "set title \"Frequency\"\n"
             "set ylabel 'ocurrences'\n"
             "set xlabel 'bg (photons)'\n"
             "set origin 0.0, 0.0\n"
             "set size 0.5, 1.0\n"
             "set grid ytics\n"
             "set format y '%%.0s%%c'\n"
             "set style fill solid 1 border rgb 'black'\n"
             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n", func_data->output_file_path);
  gnuprint( f, func_data->data, func_data->size );

  float mean = cal_mean( func_data->size, func_data->data );

  fprintf(f, "set title \"Background: Mean = %g\n"
             "set ylabel 'bg (photons)'\n"
             "set xlabel 'emissions'\n"
             "set origin 0.5, 0.0\n"
             "set size 0.5, 1.0\n"
             "set format x '%%.0s%%c'\n"
             "plot '-' with line lc rgb 'blue'\n", mean );
  gnuprint( f, func_data->data, func_data->size );

  fprintf(f, "unset multiplot\n");
  fclose(f);

  pthread_exit(NULL);

  return (void*)plotFrequencyEmissionsBg; // deal with compiler warnings
}

void *plotFrequencyEmissionsPho( void *data ) {

  struct fe_data *func_data = (struct fe_data *) data;

  FILE *f = popen("gnuplot", "w");

  fprintf(f, "set terminal png small\n" 
             //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
             "set output '%s'\n"
             "unset key\n"
             "set multiplot layout 1,2\n"
             "binwidth=50\n" 
             "bin(x,width)=width*floor(x/width)\n" 
             "set title \"Frequency\"\n"
             "set ylabel 'ocurrences'\n"
             "set xlabel 'bg (photons)'\n"
             "set origin 0.0, 0.0\n"
             "set size 0.5, 1.0\n"
             "set grid ytics\n"
             "set format y '%%.0s%%c'\n"
             "set style fill solid 1 border rgb 'black'\n"
             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n", func_data->output_file_path);

  gnuprint( f, func_data->data, func_data->size );

  float mean = cal_mean( func_data->size, func_data->data );

  fprintf(f, "set title \"Background: Mean = %g\n"
             "set ylabel 'bg (photons)'\n"
             "set xlabel 'emissions'\n"
             "set origin 0.5, 0.0\n"
             "set size 0.5, 1.0\n"
             "set format x '%%.0s%%c'\n"
             "plot '-' with line lc rgb 'blue'\n", mean );

  gnuprint( f, func_data->data, func_data->size );

  fprintf(f, "unset multiplot\n");
  fclose(f);

  return (void*)plotFrequencyEmissionsPho; // deal with compiler warnings
}

void *plotImage( void *data  ) {

  struct im_data *func_data = (struct im_data *) data;

  FILE *f = popen("gnuplot", "w");

  fprintf(f, "set terminal png small\n" 
             //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
             "set output '%s'\n"
             "unset key\n"
             "set multiplot layout 2,2\n"
             "set title 'bg'\n"
             "set palette rgbformulae 22,13,-31\n"
             "plot '-' with points palette pointtype 1\n", func_data->output_file_path);
  gnuprint( f, func_data->xtot, func_data->ytot, func_data->bgtot, func_data->size );

  fprintf(f, "set title 'photons'\n" 
             "plot '-'  u 1:($1<400 ? ($2<400?$2:1/0) : 1/0):($3<400 ? $3:1/0) with points palette pointtype 1\n");
  gnuprint( f, func_data->xtot, func_data->ytot, func_data->photot, func_data->size );

  fprintf(f, "set title 'uncertainty'\n" 
             "plot '-' with points palette pointtype 1\n");
  gnuprint( f, func_data->xtot, func_data->ytot, func_data->crlbytot, func_data->size );

  fprintf(f, "set title 'llr'\n" 
             "plot '-' with points palette pointtype 1\n");
  gnuprint( f, func_data->xtot, func_data->ytot, func_data->lltot, func_data->size );

  fprintf(f, "unset multiplot\n");
  fclose(f);

  return (void*)plotImage; // deal with compiler warnings
}

void *plotUncer( void *data ) {

  struct un_data *func_data = (struct un_data *) data;

  FILE *f = popen("gnuplot", "w");

  fprintf(f, "set terminal png small\n" 
             //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
             "set output '%s'\n"
             "unset key\n"
             "set multiplot layout 2,2\n"

             "set xrange [0:0.5]\n"
             "binwidth=0.05\n"
             "bin(x,width)=width*floor(x/width)\n"

             "set title 'uncerx Frequency'\n"
             "set ylabel 'ocurrences'\n"
             "set xlabel 'uncertainty x (photons)'\n"

             "set grid ytics\n"
             "set format y '%%.0s%%c'\n"

             "set style fill solid 1 border rgb 'black'\n"
             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n", func_data->output_file_path);
  gnuprint( f, func_data->crlbxtot, func_data->size );

  ////////////////////////////
  float xmean = cal_mean( func_data->size, func_data->crlbxtot );

  fprintf(f, "set autoscale\n"
             "set title \"uncerx: Mean = %g\n"
             "set ylabel 'uncerx (photons)'\n"
             "set xlabel 'emissions'\n"
             "set format x '%%.0s%%c'\n"
             "plot '-' with line lc rgb 'blue'\n", xmean );
  gnuprint( f, func_data->crlbxtot, func_data->size );

  ////////////////////////////
  fprintf(f, "set xrange [0:0.5]\n"
             "binwidth=0.05\n"
             "bin(x,width)=width*floor(x/width)\n"
             "set title 'uncery Frequency'\n"
             "set ylabel 'ocurrences'\n"
             "set xlabel 'uncertainty y (photons)'\n"
             "set grid ytics\n"
             "set format x '%%g'\n"
             "set format y '%%.0s%%c'\n"
             "set style fill solid 1 border rgb 'black'\n"
             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n");
  gnuprint( f, func_data->crlbytot, func_data->size );

  ///////////////////////////////

  float ymean = cal_mean( func_data->size, func_data->crlbytot );

  fprintf(f, "set autoscale\n"
             "set title \"uncery: Mean = %g\n"
             "set ylabel 'uncery (photons)'\n"
             "set xlabel 'emissions'\n"
             "set format x '%%.0s%%c'\n"
             "plot '-' with line lc rgb 'blue'\n", ymean );
  gnuprint( f, func_data->crlbytot, func_data->size );

  fprintf(f, "unset multiplot\n");
  fclose(f);

  return (void*)plotUncer; // deal with compiler warnings
}

//void plot( const char* output_path, SRscmos *sr ) {
//  // BG
//  FILE *f = popen("gnuplot", "w");
//
//  fprintf(f, "set terminal png small\n" 
//              //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
//             "set output '%s'\n"
//             "unset key\n"
//             "set multiplot layout 1,2\n"
//             "binwidth=1\n" 
//             "bin(x,width)=width*floor(x/width)\n" 
//             "set title \"Frequency\"\n"
//             "set ylabel 'ocurrences'\n"
//             "set xlabel 'bg (photons)'\n"
//             "set origin 0.0, 0.0\n"
//             "set size 0.5, 1.0\n"
//             "set grid ytics\n"
//             "set format y '%%.0s%%c'\n"
//             "set style fill solid 1 border rgb 'black'\n"
//             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n", "output/bg.png");
//  gnuprint( f, sr->bgtot, sr->total_size );
//
//  float mean = cal_mean( sr->total_size, sr->bgtot );
//
//  fprintf(f, "set title \"Background: Mean = %g\n"
//             "set ylabel 'bg (photons)'\n"
//             "set xlabel 'emissions'\n"
//             "set origin 0.5, 0.0\n"
//             "set size 0.5, 1.0\n"
//             "set format x '%%.0s%%c'\n"
//             "plot '-' with line lc rgb 'blue'\n", mean );
//  gnuprint( f, sr->bgtot, sr->total_size );
//
//  fprintf(f, "unset multiplot\n");
//  fclose(f);
//
//  // photons
//  f = popen("gnuplot", "w");
//
//  fprintf(f, "set terminal png small\n" 
//             //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
//             "set output '%s'\n"
//             "unset key\n"
//             "set multiplot layout 1,2\n"
//             "binwidth=50\n" 
//             "bin(x,width)=width*floor(x/width)\n" 
//             "set title \"Frequency\"\n"
//             "set ylabel 'ocurrences'\n"
//             "set xlabel 'bg (photons)'\n"
//             "set origin 0.0, 0.0\n"
//             "set size 0.5, 1.0\n"
//             "set grid ytics\n"
//             "set format y '%%.0s%%c'\n"
//             "set style fill solid 1 border rgb 'black'\n"
//             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n", "output/pho.png");
//
//  gnuprint( f, sr->photot, sr->total_size );
//
//  mean = cal_mean( sr->total_size, sr->photot );
//
//  fprintf(f, "set title \"Background: Mean = %g\n"
//             "set ylabel 'bg (photons)'\n"
//             "set xlabel 'emissions'\n"
//             "set origin 0.5, 0.0\n"
//             "set size 0.5, 1.0\n"
//             "set format x '%%.0s%%c'\n"
//             "plot '-' with line lc rgb 'blue'\n", mean );
//
//  gnuprint( f, sr->photot, sr->total_size );
//
//  fprintf(f, "unset multiplot\n");
//  fclose(f);
//
//  // im
//  f = popen("gnuplot", "w");
//
//  fprintf(f, "set terminal png small\n" 
//             //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
//             "set output '%s'\n"
//             "unset key\n"
//             "set multiplot layout 2,2\n"
//             "set title 'bg'\n"
//             "set palette rgbformulae 22,13,-31\n"
//             "plot '-' with points palette pointtype 1\n", "output/im.png");
//  gnuprint( f, sr->xtot, sr->ytot, sr->bgtot, sr->total_size );
//
//  fprintf(f, "set title 'photons'\n" 
//             "plot '-'  u 1:($1<400 ? ($2<400?$2:1/0) : 1/0):($3<400 ? $3:1/0) with points palette pointtype 1\n");
//  gnuprint( f, sr->xtot, sr->ytot, sr->photot, sr->total_size );
//
//  fprintf(f, "set title 'uncertainty'\n" 
//             "plot '-' with points palette pointtype 1\n");
//  gnuprint( f, sr->xtot, sr->ytot, sr->crlbytot, sr->total_size );
//
//  fprintf(f, "set title 'llr'\n" 
//             "plot '-' with points palette pointtype 1\n");
//  gnuprint( f, sr->xtot, sr->ytot, sr->lltot, sr->total_size );
//
//  fprintf(f, "unset multiplot\n");
//  fclose(f);
//
//  //crlb
//  f = popen("gnuplot", "w");
//
//  fprintf(f, "set terminal png small\n" 
//             //"set terminal pngcairo size 800,600 enhanced font 'Verdana,10'\n"
//             "set output '%s'\n"
//             "unset key\n"
//             "set multiplot layout 2,2\n"
//
//             "set xrange [0:0.5]\n"
//             "binwidth=0.05\n"
//             "bin(x,width)=width*floor(x/width)\n"
//
//             "set title 'uncerx Frequency'\n"
//             "set ylabel 'ocurrences'\n"
//             "set xlabel 'uncertainty x (photons)'\n"
//
//             "set grid ytics\n"
//             "set format y '%%.0s%%c'\n"
//
//             "set style fill solid 1 border rgb 'black'\n"
//             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n", "output/crlb.png");
//  gnuprint( f, sr->crlbxtot, sr->total_size );
//
//  ////////////////////////////
//  float xmean = cal_mean( sr->total_size, sr->crlbxtot );
//
//  fprintf(f, "set autoscale\n"
//             "set title \"uncerx: Mean = %g\n"
//             "set ylabel 'uncerx (photons)'\n"
//             "set xlabel 'emissions'\n"
//             "set format x '%%.0s%%c'\n"
//             "plot '-' with line lc rgb 'blue'\n", xmean );
//  gnuprint( f, sr->crlbxtot, sr->total_size );
//
//  ////////////////////////////
//  fprintf(f, "set xrange [0:0.5]\n"
//             "binwidth=0.05\n"
//             "bin(x,width)=width*floor(x/width)\n"
//             "set title 'uncery Frequency'\n"
//             "set ylabel 'ocurrences'\n"
//             "set xlabel 'uncertainty y (photons)'\n"
//             "set grid ytics\n"
//             "set format x '%%g'\n"
//             "set format y '%%.0s%%c'\n"
//             "set style fill solid 1 border rgb 'black'\n"
//             "plot '-' using (bin($1,binwidth)):(1.0) smooth freq with boxes lc rgb 'blue'\n");
//  gnuprint( f, sr->crlbytot, sr->total_size );
//
//  ///////////////////////////////
//
//  float ymean = cal_mean( sr->total_size, sr->crlbytot );
//
//  fprintf(f, "set autoscale\n"
//             "set title \"uncery: Mean = %g\n"
//             "set ylabel 'uncery (photons)'\n"
//             "set xlabel 'emissions'\n"
//             "set format x '%%.0s%%c'\n"
//             "plot '-' with line lc rgb 'blue'\n", ymean );
//  gnuprint( f, sr->crlbytot, sr->total_size );
//
//  fprintf(f, "unset multiplot\n");
//  fclose(f);
//}

void plot( const char* output_path, SRscmos *sr ) {
  pthread_t threads[4];
  int rc = 0; // pthreads error

  std::string plot_output = output_path;

  std::string bg_path = plot_output + "bg.png";

  struct fe_data *bg_data = (struct fe_data *)malloc( sizeof(*bg_data) );
  bg_data->output_file_path = bg_path.c_str();
  //bg_data->output_file_path = "output/bg.png";
  bg_data->data = sr->bgtot;
  bg_data->size = sr->total_size;

  std::string pho_path = plot_output + "pho.png";

  struct fe_data *pho_data = (struct fe_data *)malloc( sizeof(*pho_data) );
  pho_data->output_file_path = pho_path.c_str();
  //pho_data->output_file_path = "output/pho.png";
  pho_data->data = sr->photot;
  pho_data->size = sr->total_size;

  std::string im_path = plot_output + "im.png";

  struct im_data *implot_data = (struct im_data *) malloc( sizeof(*implot_data) );
  implot_data->output_file_path = im_path.c_str();
  //implot_data->output_file_path = "output/im.png";
  implot_data->xtot = sr->xtot;
  implot_data->ytot = sr->ytot;
  implot_data->bgtot = sr->bgtot;
  implot_data->photot = sr->photot;
  implot_data->crlbytot = sr->crlbytot;
  implot_data->lltot = sr->lltot;
  implot_data->size = sr->total_size;

  std::string crlb_path = plot_output + "crlb.png";
  
  struct un_data *cer_data = (struct un_data *) malloc( sizeof(*cer_data) );
  cer_data->output_file_path = crlb_path.c_str();
  //cer_data->output_file_path = "output/crlb.png";
  cer_data->crlbxtot = sr->crlbxtot;
  cer_data->crlbytot = sr->crlbytot;
  cer_data->size = sr->total_size;
  
  rc = pthread_create( &threads[0], NULL, &plotFrequencyEmissionsBg, bg_data );
  if (rc) {
    printf("cannot create thread\n");
  }
  rc = pthread_create( &threads[1], NULL, &plotFrequencyEmissionsPho, pho_data );
  rc = pthread_create( &threads[2], NULL, &plotImage, implot_data );
  rc = pthread_create( &threads[3], NULL, &plotUncer, cer_data );

  (void) pthread_join( threads[0], NULL );
  (void) pthread_join( threads[1], NULL );
  (void) pthread_join( threads[2], NULL );
  (void) pthread_join( threads[3], NULL );

  free(bg_data);
  free(pho_data);
  free(implot_data);
  free(cer_data);

  //rc = pthread_create( &threads[0], NULL, &plotFrequencyEmissionsBg, (void *)&thread_data[0] );

  //plotFrequencyEmissionsBg( "output/bg.png", sr->bgtot, sr->total_size, "bg", "Background", 1 );
  //plotFrequencyEmissionsPho( "output/photot.png", sr->photot, sr->total_size, "photons", "Photons", 50 );
  //plotImage( "output/im.png", sr->xtot, sr->ytot, sr->bgtot, sr->photot, sr->crlbytot, sr->lltot, sr->total_size );
  //plotUncer( "output/crlb.png", sr->crlbxtot, sr->crlbytot, sr->total_size ); 
}
