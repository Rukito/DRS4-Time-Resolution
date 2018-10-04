#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <string.h>

#include <TFile.h>
#include <TTree.h>
#include <sstream>
#include <TObject.h>

typedef struct{
	char	tag[3];
	char	version;
} FHEADER;

typedef struct {
   char           time_header[4];
} THEADER;

typedef struct {
   char           bn[2];
   unsigned short board_serial_number;
} BHEADER;

typedef struct {
   char           event_header[4];
   unsigned int   event_serial_number;
   unsigned short year;
   unsigned short month;
   unsigned short day;
   unsigned short hour;
   unsigned short minute;
   unsigned short second;
   unsigned short millisecond;
   unsigned short range;
} EHEADER;

typedef struct {
   char           tc[2];
   unsigned short trigger_cell;
} TCHEADER;

typedef struct {
   char           c[1];
   char           cn[3];
} CHEADER;

/*-----------------------------------------------------------------------------*/

int main(int argc, const char * argv[])
{
   FHEADER  fh;
   THEADER  th;
   BHEADER  bh;
   EHEADER  eh;
   TCHEADER tch;
   CHEADER  ch;
   
   unsigned int scaler;
   unsigned short voltage[1024];
   double waveform[16][4][1024], time[16][4][1024];
   double_t waveform1[1024], waveform2[1024], waveform3[1024], waveform4[1024];
   double_t time1[1024], time2[1024], time3[1024], time4[1024];
   float bin_width[16][4][1024];
   int i, j, b, chn, n, chn_index, n_boards;
   double t1, t2, dt, t3, t4, dt_;

   int ndt;
   double threshold, sumdt, sumdt2, sumdt_, sumdt_2;
   std::string root_dir, data_dir;
   
   if (argc < 2){
      printf("Usage: read_binary <path> <filename>");
      return 0;
   }
   //Create file with trees of waveforms//////////////////////////////////////////////////
   ///////////////////////////////////////////////////////////////////////////////////////

   root_dir = "../root_files/";
   for(int i=8; i<strlen(argv[1]);i++) root_dir += argv[1][i];
   root_dir += ".root";
   const char* root_fname = root_dir.c_str();
   data_dir = argv[1];
   const char* filename = argv[1];

   TFile rootf(root_fname, "recreate");
   TTree tree1("Waveforms", "Waveforms Tree");

   tree1.Branch("waveform1", &waveform1, "waveform1[1024]/D");
   tree1.Branch("waveform2", &waveform2, "waveform2[1024]/D");
   tree1.Branch("waveform3", &waveform3, "waveform3[1024]/D");
   tree1.Branch("waveform4", &waveform4, "waveform4[1024]/D");
   tree1.Branch("time1", &time1, "time1[1024]/D");
   tree1.Branch("time2", &time2, "time2[1024]/D");
   tree1.Branch("time3", &time3, "time3[1024]/D");
   tree1.Branch("time4", &time4, "time4[1024]/D");

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
  
   // open the binary waveform file
   FILE *f = fopen(filename, "rb");
   if (f == NULL) {
      printf("Cannot find file \'%s\'\n", filename);
      return 0;
   }

   // read file header
   fread(&fh, sizeof(fh), 1, f);
   if (fh.tag[0] != 'D' || fh.tag[1] != 'R' || fh.tag[2] != 'S') {
      printf("Found invalid file header in file \'%s\', aborting.\n", filename);
      return 0;
   }
   
   if (fh.version != '2') {
      printf("Found invalid file version \'%c\' in file \'%s\', should be \'2\', aborting.\n", fh.version, filename);
      return 0;
   }

   // read time header
   fread(&th, sizeof(th), 1, f);
   if (memcmp(th.time_header, "TIME", 4) != 0) {
      printf("Invalid time header in file \'%s\', aborting.\n", filename);
      return 0;
   }

   for (b = 0 ; ; b++) {
      // read board header
      fread(&bh, sizeof(bh), 1, f);
      if (memcmp(bh.bn, "B#", 2) != 0) {
         // probably event header found
         fseek(f, -4, SEEK_CUR);
         break;
      }
      
      printf("Found data for board #%d\n", bh.board_serial_number);

      // read time bin widths
      memset(bin_width[b], sizeof(bin_width[0]), 0);
      for (chn=0 ; chn<5 ; chn++) {
         fread(&ch, sizeof(ch), 1, f);
         if (ch.c[0] != 'C') {
            // event header found
            fseek(f, -4, SEEK_CUR);
            break;
         }
         i = ch.cn[2] - '0' - 1;
         printf("Found timing calibration for channel #%d\n", i+1);
         fread(&bin_width[b][i][0], sizeof(float), 1024, f);
         // fix for 2048 bin mode: double channel
         if (bin_width[b][i][1023] > 10 || bin_width[b][i][1023] < 0.01) {
            for (j=0 ; j<512 ; j++)
               bin_width[b][i][j+512] = bin_width[b][i][j];
         }
      }
   }
   n_boards = b;
   
   // initialize statistics
   ndt = 0;
   sumdt = sumdt2 = 0;
  
   // loop over all events in the data file
   for (n=0 ; ; n++) {
      // read event header
      i = (int)fread(&eh, sizeof(eh), 1, f);
      if (i < 1)
         break;
      printf("Found event #%d %d %d\n", eh.event_serial_number, eh.second, eh.millisecond);

      // loop over all boards in data file
      for (b=0 ; b<n_boards ; b++) {

         // read board header
         fread(&bh, sizeof(bh), 1, f);
         if (memcmp(bh.bn, "B#", 2) != 0) {
            printf("Invalid board header in file \'%s\', aborting.\n", filename);
            return 0;
         }

         // read trigger cell
         fread(&tch, sizeof(tch), 1, f);
         if (memcmp(tch.tc, "T#", 2) != 0) {
            printf("Invalid trigger cell header in file \'%s\', aborting.\n", filename);
            return 0;
         }

         if (n_boards > 1)
            printf("Found data for board #%d\n", bh.board_serial_number);

         // reach channel data
         for (chn=0 ; chn<4 ; chn++) {

            // read channel header
            fread(&ch, sizeof(ch), 1, f);
            if (ch.c[0] != 'C') {
               // event header found
               fseek(f, -4, SEEK_CUR);
               break;
            }
            chn_index = ch.cn[2] - '0' - 1;
            fread(&scaler, sizeof(int), 1, f);
            fread(voltage, sizeof(short), 1024, f);
            for (i=0 ; i<1024 ; i++) {
               // convert data to volts
               waveform[b][chn_index][i] = (voltage[i] / 65536. + eh.range/1000.0 - 0.5);
               // calculate time for this cell
               for (j=0,time[b][chn_index][i]=0 ; j<i ; j++)
                  time[b][chn_index][i] += bin_width[b][chn_index][(j+tch.trigger_cell) % 1024];
            }
         }

         // align cell #0 of all channels
         t1 = time[b][0][(1024-tch.trigger_cell) % 1024];
         for (chn=1 ; chn<4 ; chn++) {
            t2 = time[b][chn][(1024-tch.trigger_cell) % 1024];
            dt = t1 - t2;
            for (i=0 ; i<1024 ; i++)
               time[b][chn][i] += dt;
         }

// Filling tree with data //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
	 for(int i=0; i<1024; i++){
	    waveform1[i] = waveform[0][0][i];
	    waveform2[i] = waveform[0][1][i];
            waveform3[i] = waveform[0][2][i];
            waveform4[i] = waveform[0][3][i];
	    time1[i] = time[0][0][i];
            time2[i] = time[0][1][i];
            time3[i] = time[0][2][i];
            time4[i] = time[0][3][i];

	 }
	 tree1.Fill();
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
/*

         t1 = t2 = t3 = t4 =0;
         threshold = -0.02;

         // find peak in channel 1 above threshold
         for (i=0 ; i<1022 ; i++)
            if (waveform[b][0][i] < threshold && waveform[b][0][i+1] >= threshold) {
               t1 = (threshold-waveform[b][0][i])/(waveform[b][0][i+1]-waveform[b][0][i])*(time[b][0][i+1]-time[b][0][i])+time[b][0][i];
               break;
            }

         // find peak in channel 2 above threshold
         for (i=0 ; i<1022 ; i++)
            if (waveform[b][1][i] < threshold && waveform[b][1][i+1] >= threshold) {
               t2 = (threshold-waveform[b][1][i])/(waveform[b][1][i+1]-waveform[b][1][i])*(time[b][1][i+1]-time[b][1][i])+time[b][1][i];
               break;
            }

        // find peak in channel 3 above threshold
         for (i=0 ; i<1022 ; i++)
            if (waveform[b][2][i] < threshold && waveform[b][2][i+1] >= threshold) {
               t3 = (threshold-waveform[b][2][i])/(waveform[b][2][i+1]-waveform[b][2][i])*(time[b][2][i+1]-time[b][2][i])+time[b][2][i];
               break;
            }

         // find peak in channel 4 above threshold
         for (i=0 ; i<1022 ; i++)
            if (waveform[b][3][i] < threshold && waveform[b][3][i+1] >= threshold) {
               t4 = (threshold-waveform[b][3][i])/(waveform[b][3][i+1]-waveform[b][3][i])*(time[b][3][i+1]-time[b][3][i])+time[b][3][i];
               break;
            }

         // calculate distance of peaks with statistics
         if (t1 > 0 && t2 > 0) {
            ndt++;
            dt = t2 - t1;
            sumdt += dt;
            sumdt2 += dt*dt;
            printf("Proba: %lfns", sumdt);
         }

         // calculate distance of peaks with statistics
         if (t3 > 0 && t4 > 0) {
            ndt++;
            dt_ = t3 - t4;
            sumdt_ += dt;
            sumdt_2 += dt*dt;
            printf("Proba_: %lfns", sumdt_);
         }
*/      
      }
// Saving tree //////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
}
   tree1.Write();
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
 //  }
   return 0;
}
