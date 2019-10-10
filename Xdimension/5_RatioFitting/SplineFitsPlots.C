#include <vector>
#include "TSpline.h"

void SplineFitsPlots()
{
  gStyle->SetOptStat(0);

  TFile data_file("../3_TreePlotting/out_hists_data.root","read");
  TFile mc_file("../3_TreePlotting/out_hists_mc.root","read");

  TH1D *dataHist0;
  TH1D *mcHist0;
  TH1D *dataHist1;
  TH1D *mcHist1;
  TH1D *dataHist2;
  TH1D *mcHist2;

  TCanvas *c_data0 = (TCanvas*)data_file.Get("hist_trk_x_hit_sigma_plane_0_prof");
  TCanvas *c_mc0   = (TCanvas*)mc_file.Get("hist_trk_x_hit_sigma_plane_0_prof");
  dataHist0 = (TH1D*)c_data0->GetPrimitive("histtm_trk_x_hit_sigma_plane_0");
  mcHist0 = (TH1D*)c_mc0->GetPrimitive("histtm_trk_x_hit_sigma_plane_0");
  dataHist0->Divide(mcHist0);

  TCanvas *c_data1 = (TCanvas*)data_file.Get("hist_trk_x_hit_sigma_plane_1_prof");
  TCanvas *c_mc1   = (TCanvas*)mc_file.Get("hist_trk_x_hit_sigma_plane_1_prof");
  dataHist1 = (TH1D*)c_data1->GetPrimitive("histtm_trk_x_hit_sigma_plane_1");
  mcHist1 = (TH1D*)c_mc1->GetPrimitive("histtm_trk_x_hit_sigma_plane_1");
  dataHist1->Divide(mcHist1);

  TCanvas *c_data2 = (TCanvas*)data_file.Get("hist_trk_x_hit_sigma_plane_2_prof");
  TCanvas *c_mc2   = (TCanvas*)mc_file.Get("hist_trk_x_hit_sigma_plane_2_prof");
  dataHist2 = (TH1D*)c_data2->GetPrimitive("histtm_trk_x_hit_sigma_plane_2");
  mcHist2 = (TH1D*)c_mc2->GetPrimitive("histtm_trk_x_hit_sigma_plane_2");
  dataHist2->Divide(mcHist2);

  TSpline3* xSpline0 = new TSpline3(dataHist0,"b2e2",0,0);
  TSpline3* xSpline1 = new TSpline3(dataHist1,"b2e2",0,0);
  TSpline3* xSpline2 = new TSpline3(dataHist2,"b2e2",0,0);

  std::vector< TSpline3* > var0; var0.resize(1000);
  std::vector< TSpline3* > var1; var1.resize(1000);
  std::vector< TSpline3* > var2; var2.resize(1000);

  TRandom3* rand = new TRandom3();

  for(int i = 0; i<1000; i++){
    TH1D* temp0 = (TH1D*)dataHist0->Clone(Form("t0_%d",i));
    TH1D* temp1 = (TH1D*)dataHist1->Clone(Form("t1_%d",i));
    TH1D* temp2 = (TH1D*)dataHist2->Clone(Form("t2_%d",i));

    for(int bin = 1; bin <= temp0->GetNbinsX(); bin++){
      temp0->SetBinContent(bin, rand->Gaus(dataHist0->GetBinContent(bin),dataHist0->GetBinError(bin)));
    }
    for(int bin = 1; bin <= temp1->GetNbinsX(); bin++){
      temp1->SetBinContent(bin, rand->Gaus(dataHist1->GetBinContent(bin),dataHist1->GetBinError(bin)));
    }
    for(int bin = 1; bin <= temp2->GetNbinsX(); bin++){
      temp2->SetBinContent(bin, rand->Gaus(dataHist2->GetBinContent(bin),dataHist2->GetBinError(bin)));
    }

    var0[i] = new TSpline3(temp0,"b2e2",0,0);
    var0[i]->SetLineColor(kBlue-10);
    var1[i] = new TSpline3(temp1,"b2e2",0,0);
    var1[i]->SetLineColor(kRed-10);
    var2[i] = new TSpline3(temp2,"b2e2",0,0);
    var2[i]->SetLineColor(kGray);
  }



  TCanvas* c = new TCanvas("c");

  dataHist0->SetLineColor(kBlue);
  dataHist0->SetLineWidth(2);
  xSpline0->SetLineColor(kBlue);
  xSpline0->SetLineWidth(2);
  xSpline0->SetLineStyle(2);

  dataHist1->SetLineColor(kRed);
  dataHist1->SetLineWidth(2);
  xSpline1->SetLineColor(kRed);
  xSpline1->SetLineWidth(2);
  xSpline1->SetLineStyle(2);

  dataHist2->SetLineColor(kBlack);
  dataHist2->SetLineWidth(2);
  xSpline2->SetLineColor(kBlack);
  xSpline2->SetLineWidth(2);
  xSpline2->SetLineStyle(2);

  //dataHist0->GetYaxis()->SetTitle("Hit Charge Data/MC Ratio");
  dataHist0->GetYaxis()->SetTitle("Hit Width Data/MC Ratio");
  dataHist0->GetYaxis()->SetTitleOffset(1.2);
  dataHist0->GetXaxis()->SetTitle("Track x Position [cm]");
  dataHist0->Draw();
  /*  
  for(int i = 0; i<1000; i++){
    var0[i]->Draw("C same");
    var1[i]->Draw("C same");
    var2[i]->Draw("C same");
  }
  */
  dataHist0->Draw("same");
  dataHist1->Draw("same");
  dataHist2->Draw("same");
  xSpline0->Draw("C same");
  xSpline1->Draw("C same");
  xSpline2->Draw("C same");


  gStyle->SetOptStat(0);
  TLegend* leg = new TLegend(0.4,0.6,0.8,.8);  
  leg->AddEntry(dataHist0,"Plane 0", "l");
  leg->AddEntry(dataHist1,"Plane 1", "l");
  leg->AddEntry(dataHist2,"Plane 2", "l");
  leg->Draw("same");

  TFile* out = new TFile("out_splines.root","update");
  xSpline0->Write("trk_x_hit_Q_spline_Plane0_CV");
  xSpline1->Write("trk_x_hit_Q_spline_Plane1_CV");
  xSpline2->Write("trk_x_hit_Q_spline_Plane2_CV"); 
  
  /*
  TDirectory *dir = dir->mkdir("vars");
  dir->cd();     
  
  for(int i = 0; i < 1000; i++){
    var0[i]->Write(Form("trk_x_hit_Q_spline_Plane0_var_%d",i));
    var1[i]->Write(Form("trk_x_hit_Q_spline_Plane1_var_%d",i));
    var2[i]->Write(Form("trk_x_hit_Q_spline_Plane2_var_%d",i));
  }
  */
}
