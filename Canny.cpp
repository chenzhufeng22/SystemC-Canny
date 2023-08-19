/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "systemc.h"

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEONAME "Engineering"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */
#define SET_STACK_SIZE set_stack_size(128*1024*1024);
/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */
#define WINSIZE 21

typedef struct Tempim_s
{
	float img[SIZE];
	Tempim_s(void) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = 0;
		}
	}
	Tempim_s& operator=(const Tempim_s& copy) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = copy.img[i];
		}
		return *this;
	}
	operator float*() {
		return img;
	}
	float& operator[](const int index) {
		return img[index];
	}
}TEMPIM;

typedef struct Kernel_s
{
	float img[WINSIZE];
	Kernel_s(void) {
		for (int i = 0; i < WINSIZE; i++) {
			img[i] = 0;
		}
	}
	Kernel_s& operator=(const Kernel_s& copy) {
		for (int i = 0; i < WINSIZE; i++) {
			img[i] = copy.img[i];
		}
		return *this;
	}
	operator float*() {
		return img;
	}
	float& operator[](const int index) {
		return img[index];
	}
}KERNEL;

typedef struct Image_s
{
	unsigned char img[SIZE];

	Image_s(void)
	{
	   for (int i=0; i<SIZE; i++)
	   {
	      img[i] = 0;
	   }
	}

	Image_s& operator=(const Image_s& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   {
	      img[i] = copy.img[i];
	   }
	   return *this;
	}

	operator unsigned char*()
	{
	   return img;
	}

	unsigned char& operator[](const int index)
	{
	   return img[index];
	}
} IMAGE;

typedef struct SImage_s {
	short int img[SIZE];

	SImage_s(void) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = 0;
		}
	}

	SImage_s& operator=(const SImage_s& copy) {
		for (int i = 0; i < SIZE; i++) {
			img[i] = copy.img[i];
		}
		return *this;
	}

	operator short int*() {
		return img;
	}

	short int& operator[](const int index) {
		return img[index];
	}
} SIMAGE;

SC_MODULE(Stimulus){
	IMAGE imageout;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo_out<sc_time> timeOut;
	sc_time startTime;
	/******************************************************************************
	* Function: read_pgm_image
	* Purpose: This function reads in an image in PGM format. The image can be
	* read in from either a file or from standard input. The image is only read
	* from standard input when infilename = NULL. Because the PGM format includes
	* the number of columns and the number of rows in the image, these are read
	* from the file. Memory to store the image is allocated OUTSIDE this function.
	* The found image size is checked against the expected rows and cols.
	* All comments in the header are discarded in the process of reading the
	* image. Upon failure, this function returns 0, upon sucess it returns 1.
	******************************************************************************/
	int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
	{
	   FILE *fp;
	   char buf[71];
	   int r, c;

	   /***************************************************************************
	   * Open the input image file for reading if a filename was given. If no
	   * filename was provided, set fp to read from standard input.
	   ***************************************************************************/
	   if(infilename == NULL) fp = stdin;
	   else{
	      if((fp = fopen(infilename, "r")) == NULL){
	         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
	            infilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Verify that the image is in PGM format, read in the number of columns
	   * and rows in the image and scan past all of the header information.
	   ***************************************************************************/
	   fgets(buf, 70, fp);
	   if(strncmp(buf,"P5",2) != 0){
	      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
	   sscanf(buf, "%d %d", &c, &r);
	   if(c != cols || r != rows){
	      fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

	   /***************************************************************************
	   * Read the image from the file.
	   ***************************************************************************/
	   if((unsigned)rows != fread(image, cols, rows, fp)){
	      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }

	   if(fp != stdin) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   int i=0, n=0;
	   char infilename[40];

	   for(i=0; i<IMG_NUM; i++)
	   {
	      n = i % AVAIL_IMG;
	      sprintf(infilename, IMG_IN, n+1);

	      /****************************************************************************
	      * Read in the image.
	      ****************************************************************************/
	      if(VERBOSE) printf("Reading the image %s.\n", infilename);
	      if(read_pgm_image(infilename, imageout, ROWS, COLS) == 0){
	         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	         exit(1);
	      }
	      ImgOut.write(imageout);
				startTime = sc_time_stamp();
				timeOut.write(startTime);
				printf("%7d ms: Stimulus sent frame %02d.\n", (int)startTime.to_seconds()*1000, i+1);
	   }
	}

	SC_CTOR(Stimulus)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Monitor){
	IMAGE imagein;
	sc_fifo_in<IMAGE>  ImgIn;
	sc_fifo_in<sc_time> timeIn;
	sc_time startTime;

	/******************************************************************************
	* Function: write_pgm_image
	* Purpose: This function writes an image in PGM format. The file is either
	* written to the file specified by outfilename or to standard output if
	* outfilename = NULL. A comment can be written to the header if coment != NULL.
	******************************************************************************/
	int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
	    int cols, const char *comment, int maxval)
	{
	   FILE *fp;

	   /***************************************************************************
	   * Open the output image file for writing if a filename was given. If no
	   * filename was provided, set fp to write to standard output.
	   ***************************************************************************/
	   if(outfilename == NULL) fp = stdout;
	   else{
	      if((fp = fopen(outfilename, "w")) == NULL){
	         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
	            outfilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Write the header information to the PGM file.
	   ***************************************************************************/
	   fprintf(fp, "P5\n%d %d\n", cols, rows);
	   if(comment != NULL)
	      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
	   fprintf(fp, "%d\n", maxval);

	   /***************************************************************************
	   * Write the image data to the file.
	   ***************************************************************************/
	   if((unsigned)rows != fwrite(image, cols, rows, fp)){
	      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	      if(fp != stdout) fclose(fp);
	      return(0);
	   }

	   if(fp != stdout) fclose(fp);
	   return(1);
	}

	void main(void)
	{
		 sc_time now = SC_ZERO_TIME, last = SC_ZERO_TIME;
	   char outfilename[128];    /* Name of the output "edge" image */
	   int i, n;

	   for(i=0; i<IMG_NUM; i++)
	   {
	      ImgIn.read(imagein);
				timeIn.read(startTime);

	      /****************************************************************************
	      * Write out the edge image to a file.
	      ****************************************************************************/
	      n = i % AVAIL_IMG;
	      sprintf(outfilename, IMG_OUT, n+1);
	      if(VERBOSE) printf("Writing the edge image in the file %s.\n", outfilename);
	      if(write_pgm_image(outfilename, imagein, ROWS, COLS,"", 255) == 0){
	         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	         exit(1);
	      }
				printf("%7d ms: Monitor received frame %02d with %7d ms delay.\n",
				(int)(sc_time_stamp().to_seconds()*1000), i+1, (int)((sc_time_stamp()-startTime).to_seconds()*1000));
				last = now;
				now = sc_time_stamp();
				if(i > 0){
					printf("%7d ms: %.3f seconds after previous frame %02d, %.3f FPS.\n",
					(int)(sc_time_stamp().to_seconds()*1000), (float)(now-last).to_seconds(), i, 1/(now-last).to_seconds());
				}
	   }
	   if(VERBOSE) printf("Monitor exits simulation.\n");
	      sc_stop();	// done testing, quit the simulation
	}

	SC_CTOR(Monitor)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataIn)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}

	SC_CTOR(DataIn)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataOut)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}

	SC_CTOR(DataOut)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(RECEIVER) {
	sc_fifo_in<IMAGE> imgInP;
	sc_fifo_out<IMAGE> imgOutP;
	IMAGE img;

	void main() {
		while (1) {
			imgInP.read(img);
			imgOutP.write(img);
		}
	}

	SC_CTOR(RECEIVER) {
		SC_THREAD(main);
		SET_STACK_SIZE;
	}
};


	/*******************************************************************************
	* PROCEDURE: make_gaussian_kernel
	* PURPOSE: Create a one dimensional gaussian kernel.
	*******************************************************************************/

	SC_MODULE(GAUSSIAN_KERNEL) {
		sc_fifo_out<KERNEL> kerOutP;
		sc_fifo_out<int> wdsOutP;
		KERNEL kernel;
		int wdsize;

		void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
		{
			int i, center;
			float x, fx, sum=0.0;

			*windowsize = 1 + 2 * ceil(2.5 * sigma);
			center = (*windowsize) / 2;

			if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

			for(i=0;i<(*windowsize);i++){
				x = (float)(i - center);
				fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
				kernel[i] = fx;
				sum += fx;
			}

			for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

			if(VERBOSE){
				printf("The filter coefficients are:\n");
				for(i=0;i<(*windowsize);i++)
				printf("kernel[%d] = %f\n", i, kernel[i]);
			}
		}

		void main(){
			while (1) {
				wait(0, SC_MS);
				make_gaussian_kernel(SIGMA, kernel, &wdsize);
				wdsOutP.write(wdsize);
				kerOutP.write(kernel);
			}
		}

		SC_CTOR(GAUSSIAN_KERNEL) {
			SC_THREAD(main);
			SET_STACK_SIZE;
		}
	};


	/*******************************************************************************
	* PROCEDURE: gaussian_smooth
	* PURPOSE: Blur an image with a gaussian filter.
	*******************************************************************************/
	SC_MODULE(BLURX) {
		sc_fifo_in<IMAGE> imgInP;

		sc_fifo_in<KERNEL> kerInP;
		sc_fifo_out<KERNEL> kerOutP;

		sc_fifo_in<int> wdsInP;
		sc_fifo_out<int> wdsOutP;

		sc_fifo_out<TEMPIM> tempimOutP;

		sc_event e1, e2, e3, e4, data_received;

		IMAGE img;
		KERNEL kernel;
		int wdsize;
		TEMPIM tempim;

		void BlurXT1(){
			while(1){
				wait(data_received);
				int r, c, rr, cc,     /* Counter variables. */
				windowsize,        /* Dimension of the gaussian kernel. */
				center = wdsize /2;            /* Half of the windowsize. */
				float dot, sum;
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the x - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
				for(r = (ROWS/4)*0; r <= (rows/4)*1-1; r++){
					for(c=0;c<cols;c++){
						dot = 0.0;
						sum = 0.0;
						for(cc=(-center);cc<=center;cc++){
							if(((c+cc) >= 0) && ((c+cc) < cols)){
								dot += (float)img[r*cols+(c+cc)] * kernel[center+cc];
								sum += kernel[center+cc];
							}
						}
						tempim[r*cols+c] = dot/sum;
					}
				}
				wait(1710/4, SC_MS);
				e1.notify(SC_ZERO_TIME);
			}
		}

		void BlurXT2(){
			while(1){
				wait(data_received);
				int r, c, rr, cc,     /* Counter variables. */
				windowsize,        /* Dimension of the gaussian kernel. */
				center = wdsize /2;            /* Half of the windowsize. */
				float dot, sum;
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the x - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
				for(r = (ROWS/4)*1; r <= (rows/4)*2-1; r++){
					for(c=0;c<cols;c++){
						dot = 0.0;
						sum = 0.0;
						for(cc=(-center);cc<=center;cc++){
							if(((c+cc) >= 0) && ((c+cc) < cols)){
								dot += (float)img[r*cols+(c+cc)] * kernel[center+cc];
								sum += kernel[center+cc];
							}
						}
						tempim[r*cols+c] = dot/sum;
					}
				}
				wait(1710/4, SC_MS);
				e2.notify(SC_ZERO_TIME);
			}
		}

		void BlurXT3(){
			while(1){
				wait(data_received);
				int r, c, rr, cc,     /* Counter variables. */
				windowsize,        /* Dimension of the gaussian kernel. */
				center = wdsize /2;            /* Half of the windowsize. */
				float dot, sum;
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the x - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
				for(r = (ROWS/4)*2; r <= (rows/4)*3-1; r++){
					for(c=0;c<cols;c++){
						dot = 0.0;
						sum = 0.0;
						for(cc=(-center);cc<=center;cc++){
							if(((c+cc) >= 0) && ((c+cc) < cols)){
								dot += (float)img[r*cols+(c+cc)] * kernel[center+cc];
								sum += kernel[center+cc];
							}
						}
						tempim[r*cols+c] = dot/sum;
					}
				}
				wait(1710/4, SC_MS);
				e3.notify(SC_ZERO_TIME);
			}
		}

		void BlurXT4(){
			while(1){
				wait(data_received);
				int r, c, rr, cc,     /* Counter variables. */
				windowsize,        /* Dimension of the gaussian kernel. */
				center = wdsize /2;            /* Half of the windowsize. */
				float dot, sum;
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the x - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
				for(r = (ROWS/4)*3; r <= (rows/4)*4-1; r++){
					for(c=0;c<cols;c++){
						dot = 0.0;
						sum = 0.0;
						for(cc=(-center);cc<=center;cc++){
							if(((c+cc) >= 0) && ((c+cc) < cols)){
								dot += (float)img[r*cols+(c+cc)] * kernel[center+cc];
								sum += kernel[center+cc];
							}
						}
						tempim[r*cols+c] = dot/sum;
					}
				}
				wait(1710/4, SC_MS);
				e4.notify(SC_ZERO_TIME);
			}
		}

		void main(){
			while (1) {
				imgInP.read(img);
				kerInP.read(kernel);
				wdsInP.read(wdsize);
				data_received.notify(SC_ZERO_TIME);
				wait(e1 & e2 & e3 & e4);
				tempimOutP.write(tempim);
				wdsOutP.write(wdsize);
				kerOutP.write(kernel);
			}
		}

		SC_CTOR(BLURX) {
			SC_THREAD(main);
			SET_STACK_SIZE;
			SC_THREAD(BlurXT1);
			SET_STACK_SIZE;
			SC_THREAD(BlurXT2);
			SET_STACK_SIZE;
			SC_THREAD(BlurXT3);
			SET_STACK_SIZE;
			SC_THREAD(BlurXT4);
			SET_STACK_SIZE;
		}
	};

	SC_MODULE(BLURY) {
		sc_fifo_in<TEMPIM> tempimInP;
		sc_fifo_in<KERNEL> kerInP;
		sc_fifo_in<int> wdsInP;

		sc_fifo_out<SIMAGE> sdimOutP;

		sc_event e1, e2, e3, e4, data_received;
		TEMPIM tempim;
		KERNEL kernel;
		int wsize;
		SIMAGE smoothedim;

		void BlurYT1(){
			while(1){
				wait(data_received);
				int r, c, rr, /* Counter variables. */
				center = wsize / 2;
				float dot, sum;            /* Dot product summing variable. */
				int rows = ROWS, cols = COLS;

				/****************************************************************************
				* Blur in the y - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
				for(c = (COLS/4)*0; c <= (COLS/4)*1-1; c++){
					for(r=0;r<rows;r++){
						sum = 0.0;
						dot = 0.0;
						for(rr=(-center);rr<=center;rr++){
							if(((r+rr) >= 0) && ((r+rr) < rows)){
								dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
								sum += kernel[center+rr];
							}
						}
						smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
					}
				}
				wait(1820/4, SC_MS);
				e1.notify(SC_ZERO_TIME);
			}
		}

		void BlurYT2(){
			while(1){
				wait(data_received);
				int r, c, rr, /* Counter variables. */
				center = wsize / 2;
				float dot, sum;            /* Dot product summing variable. */
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the y - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
				for(c = (COLS/4)*1; c <= (COLS/4)*2-1; c++){
					for(r=0;r<rows;r++){
						sum = 0.0;
						dot = 0.0;
						for(rr=(-center);rr<=center;rr++){
							if(((r+rr) >= 0) && ((r+rr) < rows)){
								dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
								sum += kernel[center+rr];
							}
						}
						smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
					}
				}
				wait(1820/4, SC_MS);
				e2.notify(SC_ZERO_TIME);
			}
		}

		void BlurYT3(){
			while(1){
				wait(data_received);
				int r, c, rr, /* Counter variables. */
				center = wsize / 2;
				float dot, sum;            /* Dot product summing variable. */
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the y - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
				for(c = (COLS/4)*2; c <= (COLS/4)*3-1; c++){
					for(r=0;r<rows;r++){
						sum = 0.0;
						dot = 0.0;
						for(rr=(-center);rr<=center;rr++){
							if(((r+rr) >= 0) && ((r+rr) < rows)){
								dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
								sum += kernel[center+rr];
							}
						}
						smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
					}
				}
				wait(1820/4, SC_MS);
				e3.notify(SC_ZERO_TIME);
			}
		}

		void BlurYT4(){
			while(1){
				wait(data_received);
				int r, c, rr, /* Counter variables. */
				center = wsize / 2;
				float dot, sum;            /* Dot product summing variable. */
				int rows = ROWS, cols = COLS;
				/****************************************************************************
				* Blur in the y - direction.
				****************************************************************************/
				if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
				for(c = (COLS/4)*3; c <= (COLS/4)*4-1; c++){
					for(r=0;r<rows;r++){
						sum = 0.0;
						dot = 0.0;
						for(rr=(-center);rr<=center;rr++){
							if(((r+rr) >= 0) && ((r+rr) < rows)){
								dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
								sum += kernel[center+rr];
							}
						}
						smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
					}
				}
				wait(1820/4, SC_MS);
				e4.notify(SC_ZERO_TIME);
			}
		}



		void main(){
			while (1) {
				tempimInP.read(tempim);
				kerInP.read(kernel);
				wdsInP.read(wsize);
				data_received.notify(SC_ZERO_TIME);
				wait(e1 & e2 & e3 & e4);
				sdimOutP.write(smoothedim);
			}
		}

		SC_CTOR(BLURY) {
			SC_THREAD(main);
			SET_STACK_SIZE;
			SC_THREAD(BlurYT1);
			SET_STACK_SIZE;
			SC_THREAD(BlurYT2);
			SET_STACK_SIZE;
			SC_THREAD(BlurYT3);
			SET_STACK_SIZE;
			SC_THREAD(BlurYT4);
			SET_STACK_SIZE;
		}
	};


	SC_MODULE(GAUSSIAN_SMOOTH) {
	sc_fifo_in<IMAGE> imgInP;
	sc_fifo_out<SIMAGE> sdimOutP;

	sc_fifo<IMAGE> imgCh;
	sc_fifo<KERNEL> kerChX, kerChY;
	sc_fifo<int> wdsizeChX, wdsizeChY;
	sc_fifo<TEMPIM> tepimCh;

	BLURX blurXMod;
	BLURY blurYMod;
	RECEIVER recMod;
	GAUSSIAN_KERNEL kerMod;


	void before_end_of_elaboration() {
		recMod.imgInP.bind(imgInP);
		recMod.imgOutP.bind(imgCh);

		kerMod.kerOutP.bind(kerChX);
		kerMod.wdsOutP.bind(wdsizeChX);

		blurXMod.imgInP.bind(imgCh);
		blurXMod.wdsInP.bind(wdsizeChX);
		blurXMod.kerInP.bind(kerChX);
		blurXMod.kerOutP.bind(kerChY);
		blurXMod.tempimOutP.bind(tepimCh);
		blurXMod.wdsOutP.bind(wdsizeChY);

		blurYMod.kerInP.bind(kerChY);
		blurYMod.tempimInP.bind(tepimCh);
		blurYMod.sdimOutP.bind(sdimOutP);
		blurYMod.wdsInP.bind(wdsizeChY);
	}


	SC_CTOR(GAUSSIAN_SMOOTH) :
		imgCh("imgCh", 1), kerChX("kerChX", 1), kerChY("kerChY", 1),
		wdsizeChX("wdsizeChX", 1), wdsizeChY("wdsizeChY", 1),
		tepimCh("tepimCh", 1), recMod("recMod"), kerMod("kerMod"),
		blurXMod("blurXMod"), blurYMod("blurYMod") {
	}
};


	/*******************************************************************************
	* PROCEDURE: derivative_x_y
	* PURPOSE: Compute the first derivative of the image in both the x any y
	* directions. The differential filters that are used are:
	*
	*                                          -1
	*         dx =  -1 0 +1     and       dy =  0
	*                                          +1
	*
	*******************************************************************************/

	SC_MODULE(DERIVATIVE_X_Y) {
		sc_fifo_in<SIMAGE> sdimInP;
		sc_fifo_out<SIMAGE> dxOutP, dyOutP;
		SIMAGE stdim, outX, outY;

		void derivative_x_y(short int *smoothedim, int rows, int cols,
		        short int *delta_x, short int *delta_y)
		{
		   int r, c, pos;

		   /****************************************************************************
		   * Compute the x-derivative. Adjust the derivative at the borders to avoid
		   * losing pixels.
		   ****************************************************************************/
		   if(VERBOSE) printf("   Computing the X-direction derivative.\n");
		   for(r=0;r<rows;r++){
		      pos = r * cols;
		      delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
		      pos++;
		      for(c=1;c<(cols-1);c++,pos++){
		         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
		      }
		      delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
		   }

		   /****************************************************************************
		   * Compute the y-derivative. Adjust the derivative at the borders to avoid
		   * losing pixels.
		   ****************************************************************************/
		   if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
		   for(c=0;c<cols;c++){
		      pos = c;
		      delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
		      pos += cols;
		      for(r=1;r<(rows-1);r++,pos+=cols){
		         delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
		      }
		      delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
		   }
		}

		void main(){
			while(1){
				sdimInP.read(stdim);
				wait(480, SC_MS);
				derivative_x_y(stdim, ROWS, COLS, outX, outY);
				dxOutP.write(outX);
				dyOutP.write(outY);
			}
		}

		SC_CTOR(DERIVATIVE_X_Y) {
			SC_THREAD(main);
			SET_STACK_SIZE;
		}
	};



	/*******************************************************************************
	* PROCEDURE: magnitude_x_y
	* PURPOSE: Compute the magnitude of the gradient. This is the square root of
	* the sum of the squared derivative values.
	*******************************************************************************/
	SC_MODULE(MAGNITUDE_X_Y) {

		sc_fifo_in<SIMAGE> dxInP, dyInP;
		sc_fifo_out<SIMAGE> magOutP, dxOutP, dyOutP;
		SIMAGE dx, dy, mag;

		void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
			short int *magnitude)
			{
				int r, c, pos, sq1, sq2;

				for(r=0,pos=0;r<rows;r++){
					for(c=0;c<cols;c++,pos++){
						sq1 = (int)delta_x[pos] * (int)delta_x[pos];
						sq2 = (int)delta_y[pos] * (int)delta_y[pos];
						magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
					}
				}

			}

			void main(){
				while (1) {
					dxInP.read(dx);
					dyInP.read(dy);
					dxOutP.write(dx);
					dyOutP.write(dy);
					wait(1030, SC_MS);
					magnitude_x_y(dx, dy, ROWS, COLS, mag);
					magOutP.write(mag);
				}
			}

			SC_CTOR(MAGNITUDE_X_Y) {
				SC_THREAD(main);
				SET_STACK_SIZE;
			}
		};








	/*******************************************************************************
	* PROCEDURE: non_max_supp
	* PURPOSE: This routine applies non-maximal suppression to the magnitude of
	* the gradient image.
	*******************************************************************************/
	SC_MODULE(NON_MAX_SUPP) {
		sc_fifo_in<SIMAGE> magInP, dxInP, dyInP;
		sc_fifo_out<IMAGE> nmsOutP;
		sc_fifo_out<SIMAGE>magOutP;

		SIMAGE mag, dx, dy;
		IMAGE nms;

		void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
			unsigned char *result)
			{
				int rowcount, colcount,count;
				short *magrowptr,*magptr;
				short *gxrowptr,*gxptr;
				short *gyrowptr,*gyptr,z1,z2;
				short m00,gx,gy;
				float mag1,mag2,xperp,yperp;
				unsigned char *resultrowptr, *resultptr;

				/****************************************************************************
				* Zero the edges of the result image.
				****************************************************************************/
				for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
				count<ncols; resultptr++,resultrowptr++,count++){
					*resultrowptr = *resultptr = (unsigned char) 0;
				}

				for(count=0,resultptr=result,resultrowptr=result+ncols-1;
					count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
						*resultptr = *resultrowptr = (unsigned char) 0;
					}

					/****************************************************************************
					* Suppress non-maximum points.
					****************************************************************************/
					for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
						gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
						rowcount<=nrows-2;	// bug fix 3/29/17, RD
						rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
						resultrowptr+=ncols){
							for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
								resultptr=resultrowptr;colcount<=ncols-2;	// bug fix 3/29/17, RD
								colcount++,magptr++,gxptr++,gyptr++,resultptr++){
									m00 = *magptr;
									if(m00 == 0){
										*resultptr = (unsigned char) NOEDGE;
									}
									else{
										xperp = -(gx = *gxptr)/((float)m00);
										yperp = (gy = *gyptr)/((float)m00);
									}

									if(gx >= 0){
										if(gy >= 0){
											if (gx >= gy)
											{
												/* 111 */
												/* Left point */
												z1 = *(magptr - 1);
												z2 = *(magptr - ncols - 1);

												mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

												/* Right point */
												z1 = *(magptr + 1);
												z2 = *(magptr + ncols + 1);

												mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
											}
											else
											{
												/* 110 */
												/* Left point */
												z1 = *(magptr - ncols);
												z2 = *(magptr - ncols - 1);

												mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

												/* Right point */
												z1 = *(magptr + ncols);
												z2 = *(magptr + ncols + 1);

												mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp;
											}
										}
										else
										{
											if (gx >= -gy)
											{
												/* 101 */
												/* Left point */
												z1 = *(magptr - 1);
												z2 = *(magptr + ncols - 1);

												mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

												/* Right point */
												z1 = *(magptr + 1);
												z2 = *(magptr - ncols + 1);

												mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
											}
											else
											{
												/* 100 */
												/* Left point */
												z1 = *(magptr + ncols);
												z2 = *(magptr + ncols - 1);

												mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

												/* Right point */
												z1 = *(magptr - ncols);
												z2 = *(magptr - ncols + 1);

												mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp;
											}
										}
									}
									else
									{
										if ((gy = *gyptr) >= 0)
										{
											if (-gx >= gy)
											{
												/* 011 */
												/* Left point */
												z1 = *(magptr + 1);
												z2 = *(magptr - ncols + 1);

												mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

												/* Right point */
												z1 = *(magptr - 1);
												z2 = *(magptr + ncols - 1);

												mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
											}
											else
											{
												/* 010 */
												/* Left point */
												z1 = *(magptr - ncols);
												z2 = *(magptr - ncols + 1);

												mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

												/* Right point */
												z1 = *(magptr + ncols);
												z2 = *(magptr + ncols - 1);

												mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
											}
										}
										else
										{
											if (-gx > -gy)
											{
												/* 001 */
												/* Left point */
												z1 = *(magptr + 1);
												z2 = *(magptr + ncols + 1);

												mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

												/* Right point */
												z1 = *(magptr - 1);
												z2 = *(magptr - ncols - 1);

												mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
											}
											else
											{
												/* 000 */
												/* Left point */
												z1 = *(magptr + ncols);
												z2 = *(magptr + ncols + 1);

												mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

												/* Right point */
												z1 = *(magptr - ncols);
												z2 = *(magptr - ncols - 1);

												mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
											}
										}
									}

									/* Now determine if the current point is a maximum point */

									if ((mag1 > 0.0) || (mag2 > 0.0))
									{
										*resultptr = (unsigned char) NOEDGE;
									}
									else
									{
										if (mag2 == 0.0)
										*resultptr = (unsigned char) NOEDGE;
										else
										*resultptr = (unsigned char) POSSIBLE_EDGE;
									}
								}
							}
						}
						void main(){
							while (1) {
								magInP.read(mag);
								dxInP.read(dx);
								dyInP.read(dy);
								wait(830, SC_MS);
								non_max_supp(mag, dx, dy, ROWS, COLS, nms);
								nmsOutP.write(nms);
								magOutP.write(mag);
							}
						}

						SC_CTOR(NON_MAX_SUPP) {
							SC_THREAD(main);
							SET_STACK_SIZE;
						}
					};



	/*******************************************************************************
	* PROCEDURE: apply_hysteresis
	* PURPOSE: This routine finds edges that are above some high threshhold or
	* are connected to a high pixel by a path of pixels greater than a low
	* threshold.
	*******************************************************************************/
	SC_MODULE(APPLY_HYSTERESIS) {
		sc_fifo_in<SIMAGE> magInP;
		sc_fifo_in<IMAGE> nmsInP;
		sc_fifo_out<IMAGE> edgeOutP;
		SIMAGE mag;
		IMAGE nms, edge;

		void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
			float tlow, float thigh, unsigned char *edge)
			{
				int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
				short int maximum_mag;

				/****************************************************************************
				* Initialize the edge map to possible edges everywhere the non-maximal
				* suppression suggested there could be an edge except for the border. At
				* the border we say there can not be an edge because it makes the
				* follow_edges algorithm more efficient to not worry about tracking an
				* edge off the side of the image.
				****************************************************************************/
				for(r=0,pos=0;r<rows;r++){
					for(c=0;c<cols;c++,pos++){
						if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
						else edge[pos] = NOEDGE;
					}
				}

				for(r=0,pos=0;r<rows;r++,pos+=cols){
					edge[pos] = NOEDGE;
					edge[pos+cols-1] = NOEDGE;
				}
				pos = (rows-1) * cols;
				for(c=0;c<cols;c++,pos++){
					edge[c] = NOEDGE;
					edge[pos] = NOEDGE;
				}

				/****************************************************************************
				* Compute the histogram of the magnitude image. Then use the histogram to
				* compute hysteresis thresholds.
				****************************************************************************/
				for(r=0;r<32768;r++) hist[r] = 0;
				for(r=0,pos=0;r<rows;r++){
					for(c=0;c<cols;c++,pos++){
						if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
					}
				}

				/****************************************************************************
				* Compute the number of pixels that passed the nonmaximal suppression.
				****************************************************************************/
				for(r=1,numedges=0;r<32768;r++){
					if(hist[r] != 0) maximum_mag = r;
					numedges += hist[r];
				}

				highcount = (int)(numedges * thigh + 0.5);

				/****************************************************************************
				* Compute the high threshold value as the (100 * thigh) percentage point
				* in the magnitude of the gradient histogram of all the pixels that passes
				* non-maximal suppression. Then calculate the low threshold as a fraction
				* of the computed high threshold value. John Canny said in his paper
				* "A Computational Approach to Edge Detection" that "The ratio of the
				* high to low threshold in the implementation is in the range two or three
				* to one." That means that in terms of this implementation, we should
				* choose tlow ~= 0.5 or 0.33333.
				****************************************************************************/
				r = 1;
				numedges = hist[1];
				while((r<(maximum_mag-1)) && (numedges < highcount)){
					r++;
					numedges += hist[r];
				}
				highthreshold = r;
				lowthreshold = (int)(highthreshold * tlow + 0.5);

				if(VERBOSE){
					printf("The input low and high fractions of %f and %f computed to\n",
					tlow, thigh);
					printf("magnitude of the gradient threshold values of: %d %d\n",
					lowthreshold, highthreshold);
				}

				/****************************************************************************
				* This loop looks for pixels above the highthreshold to locate edges and
				* then calls follow_edges to continue the edge.
				****************************************************************************/
				for(r=0,pos=0;r<rows;r++){
					for(c=0;c<cols;c++,pos++){
						if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
							edge[pos] = EDGE;
							follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
						}
					}
				}

				/****************************************************************************
				* Set all the remaining possible edges to non-edges.
				****************************************************************************/
				for(r=0,pos=0;r<rows;r++){
					for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
				}
			}
			/*******************************************************************************
			* PROCEDURE: follow_edges
			* PURPOSE: This procedure edges is a recursive routine that traces edgs along
			* all paths whose magnitude values remain above some specifyable lower
			* threshhold.
			*******************************************************************************/
			void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
				int cols)
				{
					short *tempmagptr;
					unsigned char *tempmapptr;
					int i;
					int x[8] = {1,1,0,-1,-1,-1,0,1},
					y[8] = {0,1,1,1,0,-1,-1,-1};

					for(i=0;i<8;i++){
						tempmapptr = edgemapptr - y[i]*cols + x[i];
						tempmagptr = edgemagptr - y[i]*cols + x[i];

						if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
							*tempmapptr = (unsigned char) EDGE;
							follow_edges(tempmapptr,tempmagptr, lowval, cols);
						}
					}
				}

				void main(){
					while (1) {
						magInP.read(mag);
						nmsInP.read(nms);
						wait(670, SC_MS);
						apply_hysteresis(mag, nms, ROWS, COLS, TLOW, THIGH, edge);
						edgeOutP.write(edge);
					}
				}

				SC_CTOR(APPLY_HYSTERESIS) {
					SC_THREAD(main);
					SET_STACK_SIZE;
				}
			};

	/*******************************************************************************
	* PROCEDURE: canny
	* PURPOSE: To perform canny edge detection.
	*******************************************************************************/
	// void canny(unsigned char *image, int rows, int cols, float sigma,
	//          float tlow, float thigh, unsigned char *edge)
	// {
	//    unsigned char nms[SIZE]    /* Points that are local maximal magnitude. */
	// 			= {0};
	//    short int smoothedim[SIZE] /* The image after gaussian smoothing.      */
	// 			= {0},
	//              delta_x[SIZE]    /* The first devivative image, x-direction. */
	// 			= {0},
	//              delta_y[SIZE]    /* The first derivative image, y-direction. */
	// 			= {0},
	//              magnitude[SIZE]  /* The magnitude of the gadient image.      */
	// 			= {0};
	//
	//    /****************************************************************************
	//    * Perform gaussian smoothing on the image using the input standard
	//    * deviation.
	//    ****************************************************************************/
	//    if(VERBOSE) printf("Smoothing the image using a gaussian kernel.\n");
	//    gaussian_smooth(image, rows, cols, sigma, smoothedim);
	//
	//    /****************************************************************************
	//    * Compute the first derivative in the x and y directions.
	//    ****************************************************************************/
	//    if(VERBOSE) printf("Computing the X and Y first derivatives.\n");
	//    derivative_x_y(smoothedim, rows, cols, delta_x, delta_y);
	//
	//    /****************************************************************************
	//    * Compute the magnitude of the gradient.
	//    ****************************************************************************/
	//    if(VERBOSE) printf("Computing the magnitude of the gradient.\n");
	//    magnitude_x_y(delta_x, delta_y, rows, cols, magnitude);
	//
	//    /****************************************************************************
	//    * Perform non-maximal suppression.
	//    ****************************************************************************/
	//    if(VERBOSE) printf("Doing the non-maximal suppression.\n");
	//    non_max_supp(magnitude, delta_x, delta_y, rows, cols, nms);
	//
	//    /****************************************************************************
	//    * Use hysteresis to mark the edge pixels.
	//    ****************************************************************************/
	//    if(VERBOSE) printf("Doing hysteresis thresholding.\n");
	//    apply_hysteresis(magnitude, nms, rows, cols, tlow, thigh, edge);
	// }

	SC_MODULE(DUT) {
		sc_fifo_in<IMAGE> imgInP;
		sc_fifo_out<IMAGE> imgOutP;

		sc_fifo<SIMAGE> sdimCh;
		sc_fifo<SIMAGE> xMagCh;
		sc_fifo<SIMAGE> yMagCh;
		sc_fifo<SIMAGE> nonXCh;
		sc_fifo<SIMAGE> nonYCh;
		sc_fifo<SIMAGE> nonMagCh;
		sc_fifo<SIMAGE> appMagCh;
		sc_fifo<IMAGE> nmsCh;

		GAUSSIAN_SMOOTH g_sMod;
		DERIVATIVE_X_Y d_xyMod;
		MAGNITUDE_X_Y mag_xyMod;
		NON_MAX_SUPP n_m_sMod;
		APPLY_HYSTERESIS a_hMod;

		void before_end_of_elaboration() {
			g_sMod.imgInP.bind(imgInP);
			a_hMod.edgeOutP.bind(imgOutP);

			g_sMod.sdimOutP.bind(sdimCh);
			d_xyMod.sdimInP.bind(sdimCh);

			d_xyMod.dxOutP.bind(xMagCh);
			mag_xyMod.dxInP.bind(xMagCh);

			d_xyMod.dyOutP.bind(yMagCh);
			mag_xyMod.dyInP.bind(yMagCh);

			mag_xyMod.magOutP.bind(nonMagCh);
			n_m_sMod.magInP.bind(nonMagCh);

			mag_xyMod.dxOutP.bind(nonXCh);
			n_m_sMod.dxInP.bind(nonXCh);

			mag_xyMod.dyOutP.bind(nonYCh);
			n_m_sMod.dyInP.bind(nonYCh);

			n_m_sMod.nmsOutP.bind(nmsCh);
			a_hMod.nmsInP.bind(nmsCh);

			n_m_sMod.magOutP.bind(appMagCh);
			a_hMod.magInP.bind(appMagCh);
		}

		SC_CTOR(DUT) :
		sdimCh("sdimCh", 1), xMagCh("xMagCh", 1), yMagCh("yMagCh", 1),
		nonXCh("nonXCh", 1), nonYCh("nonYCh", 1), nonMagCh("nonMagCh", 1),
		appMagCh("appMagCh", 1), nmsCh("nmsCh", 1),
		g_sMod("g_sMod"), d_xyMod("d_xyMod"),
		mag_xyMod("mag_xyMod"), n_m_sMod("n_m_sMod"),
		a_hMod("a_hMod") {

		}
	};

SC_MODULE(Platform)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;

	DataIn din;
	DUT canny;
	DataOut dout;

	void before_end_of_elaboration(){
	   din.ImgIn.bind(ImgIn);
	   din.ImgOut.bind(q1);
	   canny.imgInP.bind(q1);
	   canny.imgOutP.bind(q2);
	   dout.ImgIn.bind(q2);
	   dout.ImgOut.bind(ImgOut);
	}

	SC_CTOR(Platform)
	:q1("q1",2)
	,q2("q2",2)
	,din("din")
	,canny("canny")
	,dout("dout")
	{}
};

SC_MODULE(Top)
{
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;
	Stimulus stimulus;
	Platform platform;
	Monitor monitor;
	sc_fifo<sc_time> timeCh;

	void before_end_of_elaboration(){
	   stimulus.ImgOut.bind(q1);
		 stimulus.timeOut.bind(timeCh);
	   platform.ImgIn.bind(q1);
	   platform.ImgOut.bind(q2);
	   monitor.ImgIn.bind(q2);
		 monitor.timeIn.bind(timeCh);
	}

	SC_CTOR(Top)
	:q1("q1",2)
	,q2("q2",2)
	,stimulus("stimulus")
	,platform("platform")
	,monitor("monitor")
	{}
};

Top top("top");

int sc_main(int argc, char* argv[])
{
	sc_start();
	return 0;
}
