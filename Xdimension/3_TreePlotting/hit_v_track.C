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
  
  std::vector<string> kin     = {"trk_x"};
  std::vector<string> kin_lab = {"Track x Position [cm]"};
  
  double xbins0[28] = {1.6, 3.0904, 4.8292, 6.8164, 9.052, 11.536, 14.5168, 17.9944, 21.9688, 26.6884, 32.4016, 39.6052, 48.5476, 60.4708, 77.362, 102.699, 135.736, 164.302, 184.671, 199.078, 210.256, 218.95, 226.154, 232.115, 237.332, 241.803, 245.529, 250};
  double xbins1[30] = {1.6, 2.842, 4.3324, 6.0712, 7.81, 9.7972, 12.0328, 14.7652, 17.9944, 21.7204, 25.9432, 31.1596, 37.3696, 45.07, 55.006, 68.4196, 88.0432, 117.106, 149.398, 173.741, 190.881, 203.549, 213.237, 221.186, 227.644, 233.109, 237.828, 241.803, 245.28, 250};
  double xbins2[18] = {1.6, 4.3324, 7.5616, 11.536, 16.7524, 23.4592, 32.4016, 44.8216, 63.2032, 93.508, 138.717, 175.728, 198.83, 214.479, 225.905, 234.599, 241.554, 250};

  

  std::vector<double*> xbins; 
  xbins.push_back(xbins0);
  xbins.push_back(xbins1);
  xbins.push_back(xbins2);
  
  double Nbin[3] = {27,29,17};

  std::vector< std::vector<double*> > kin_bins = {  {xbins[0]},{xbins[1]},{xbins[2]}};
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
				  500, hitpr_low[h], hitpr_high[h]);

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



 TFile* out = new TFile("out_hists_mc_CVhighstats.root","RECREATE");
 
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
