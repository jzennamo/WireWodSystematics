#define hit_v_track_cxx
#include "hit_v_track.h"
#include <iostream>
#include <TH2.h>
#include <TH1.h>
#include <TH1D.h>
#include <chrono>
#include <thread>
#include <vector>
#include <TStyle.h>
#include <TCanvas.h>
#include "TObject.h"
#include "TString.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TVectorD.h"
#include <string>
#include <vector>
#include <map>


void hit_v_track::Loop()
{

  std::vector< std::vector < std::vector <TH2D*> > > hist;
  std::vector< std::vector < std::vector <TH1F*> > > histtm;
  
  std::vector<string> kin     = {"trk_x","trk_ThetaXZ"};
  std::vector<string> kin_lab = {"Track x Position [cm]","Track #theta_{XZ}"};

  std::vector<int> Nbin;
  Nbin.resize(kin.size());
  
  double xbins0[18] = {1.6,4.3324,7.3132,11.0392,15.7588,21.9688,30.166,41.344,58.2352,87.0496,135.984,176.722,200.568,216.218,227.396,236.09,242.796,250};
  double xbins1[19] = {1.6,4.084,6.8164,10.0456,14.02,18.988,25.4464,34.1404,46.312,64.942,98.2276,149.15,183.677,204.791,218.702,228.886,236.835,243.045,250};
  double xbins2[11] = {1.6,6.568,13.2748,23.2108,39.1084,68.1712,130.023,185.913,214.23,231.37,250};

  std::vector<double*> xbins; 
  xbins.push_back(xbins0);
  xbins.push_back(xbins1);
  xbins.push_back(xbins2);

  Nbin[0] = {17,18,10};

  double thetaxzbins0[19]= {-3.2, -2.2144, -1.8816, -1.6512, -1.4336, -1.2544, -1.0752, -0.896, -0.6976, -0.416, 0.3776, 0.6976, 0.9024, 1.088, 1.2608, 1.44, 1.6576, 1.9136, 3.2};
  double thetaxzbins1[20]= {-3.2, -2.5536, -2.2976, -2.08, -1.888, -1.7216, -1.5168, -1.3184, -1.088, -0.8, 0.1536, 0.8576, 1.1456, 1.3568, 1.5552, 1.7536, 1.9392, 2.1376, 2.3616, 3.2};
  double thetaxzbins2[12]= {-3.2, -2.3872, -1.9136, -1.3248, -0.9472, -0.6144, 0.1088, 0.672, 1.0112, 1.4144, 2.0544, 3.2 };

  std::vector<double*> thetaxzbins; 
  thetaxzbins.push_back(thetaxzbins0);
  thetaxzbins.push_back(thetaxzbins1);
  thetaxzbins.push_back(thetaxzbins2);

  Nbin[1] = {18,19,11};


  std::vector< std::vector<double*> > kin_bins = {  {xbins[0],xbins[1],xbins[2]}, {thetaxzbins[0],thetaxzbins[1],thetaxzbins[2]}};
  std::vector<string> hitpr     =   {"hit_Q","hit_sigma"};
  std::vector<string> hitpr_lab =   {"Integrate Hit Charge [ADC*ticks]","Hit Width [ticks]"};
  std::vector<double> hitpr_low =   {      0,          0};
  std::vector<double> hitpr_high=   {   1000,        20};
  std::vector<double> hitpr_high_plot =   {   1000,        20};
  std::vector<string> proj = {"colz", "proj_x","proj_y","prof","over"};

  // Plane //hitpr //kin //kin bin // value
  std::vector< std::vector< std::vector< std::vector < std::vector<double> > > > >  values;
  

  // Let's initialize everything
  hist.resize(3);
  histtm.resize(3);
  std::cout << "size " << histtm.size() << std::endl; 
  values.resize(3);
  for(int pl = 0; pl < 3; pl++){
    hist[pl].resize(kin.size());
    histtm[pl].resize(kin.size());
    std::cout << "size " << histtm[pl].size() << std::endl; 
    values[pl].resize(kin.size());
    for(int k = 0; k < kin.size(); k++){
      histtm[pl][k].resize(hitpr.size());
      hist[pl][k].resize(hitpr.size());
      std::cout << "size " << histtm[pl][k].size() << std::endl; 
      values[pl][k].resize(hitpr.size());
      for(int h = 0; h < hitpr.size(); h++){
	hist[pl][k][h] = new TH2D(Form("hist_%s_%s_plane_%d", kin[k].c_str(), hitpr[h].c_str(),pl),
				  Form(";%s;%s",kin_lab[k].c_str(),hitpr_lab[h].c_str()),
				  Nbin[pl],kin_bins[pl][k],				  
				  100, hitpr_low[h], hitpr_high[h]);

	histtm[pl][k][h] = new TH1F(Form("histtm_%s_%s_plane_%d", kin[k].c_str(), hitpr[h].c_str(),pl),
				    Form(";%s;%s",kin_lab[k].c_str(),hitpr_lab[h].c_str()),
				    Nbin[pl],kin_bins[pl][k]);
				    //1000, hitpr_low[h], hitpr_high[h]);

	
	values[pl][k][h].resize(Nbin[pl]);

	//std::fill(values[pl][k][h].begin(), values[pl][k][h].end(), 0);

      }// loop over hitpr
    }//loop over kin
  }//loop over planes

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(max(trk_StartX,trk_EndX) < 256 && min(trk_StartX,trk_EndX) > 0.3 && abs(hit_sigma - 5.333333) > 0.0001/* &&  trk_fracMC > 0.5*/){      
	
	// trk x // hit Q
	hist[hit_plane][0][0]->Fill(trk_x,hit_Q);	  

	// trk x // hit Q
	hist[hit_plane][0][1]->Fill(trk_x,hit_sigma);	  
	
	// Plane //hitpr //kin //kin bin // value                                                                                                          
	
	for(int i = 0; i < Nbin[hit_plane]; i++){
	  
	  if( trk_x > xbins[hit_plane][i] && trk_x < xbins[hit_plane][i+1]){
	    values[hit_plane][0][0][i].push_back(hit_Q);
	    values[hit_plane][0][1][i].push_back(hit_sigma);
	  }// check we're in the correct bin
	}// iterate through the bins 
   	
      } // cuts for analysis
   }//Iterate through tree



 TFile* out = new TFile("out_hists_data.root","RECREATE");
 
 for(int pl = 0; pl < 3; pl++){
   for(int k = 0; k < kin.size(); k++){  
     for(int h = 0; h < hitpr.size(); h++){
       


	/// Plot the colz plot 
	TCanvas* canv_0 = new TCanvas(Form("hist_%s_%s_plane_%d_%s", kin[k].c_str(), hitpr[h].c_str(),pl,proj[0].c_str()));
	canv_0->cd();
	hist[pl][k][h]->GetYaxis()->SetRangeUser(hitpr_low[h],hitpr_high_plot[h]);
	hist[pl][k][h]->Draw("colz");
	canv_0->Write();
	canv_0->Close();
	/*
	// Plot the X-projection plot
	TCanvas* canv_1 = new TCanvas(Form("hist_%s_%s_plane_%d_%s", kin[k].c_str(), hitpr[h].c_str(),pl,proj[1].c_str()));
	canv_1->cd();
	hist[pl][k][h]->ProjectionX()->Draw("");
	canv_1->Write();
	canv_1->Close();

	// Plot the Y-projection plot
	TCanvas* canv_2 = new TCanvas(Form("hist_%s_%s_plane_%d_%s", kin[k].c_str(), hitpr[h].c_str(),pl,proj[2].c_str()));
	canv_2->cd();
	hist[pl][k][h]->GetXaxis()->SetRangeUser(hitpr_low[h],hitpr_high_plot[h]);
	hist[pl][k][h]->ProjectionY()->Draw("");
	canv_2->Write();
	canv_2->Close();
	*/
	std::cout <<Form("hist_%s_%s_plane_%d_%s", kin[k].c_str(), hitpr[h].c_str(),pl,proj[0].c_str()) << std::endl;

	for(int bin = 1; bin <= histtm[pl][k][h]->GetNbinsX(); bin++){
	  
	  double tm = hit_v_track::TruncMean(values[pl][k][h][bin-1]);
	  double err = tm * (1/sqrt(values[pl][k][h][bin-1].size()));
	  
	  std::cout << tm << std::endl;
	  
 	  histtm[pl][k][h]->SetBinContent(bin, tm);
	  histtm[pl][k][h]->SetBinError(bin, err);

	}


	// Plot the Y-projection plot
	TCanvas* canv_3 = new TCanvas(Form("hist_%s_%s_plane_%d_%s", kin[k].c_str(), hitpr[h].c_str(),pl,proj[3].c_str()));
	canv_3->cd();
	histtm[pl][k][h]->Draw("");
	canv_3->Write();
	canv_3->Close();
	
	// Plot the profile over the colz
	TCanvas* canv_4 = new TCanvas(Form("hist_%s_%s_plane_%d_%s", kin[k].c_str(), hitpr[h].c_str(),pl,proj[4].c_str()));
	canv_4->cd();
	//	hist[pl][k][h]->GetYaxis()->SetRangeUser(hitpr_low[h],hitpr_high_plot[h]);
	hist[pl][k][h]->Draw("colz");
	histtm[pl][k][h]->Draw("same");
	canv_4->Write();
	canv_4->Close();

     }
   }
 }


}// game over man.
