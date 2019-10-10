void PlaneWiseRatioPlots()
{
  gStyle->SetOptStat(0);

  TFile data_file("../3_TreePlotting/out_hists_mc_DL.root","read");
  TFile mc_file("../3_TreePlotting/out_hists_mc_CVhighstats.root","read");

  TH1D *dataHist0;
  TH1D *mcHist0;
  TH1D *dataHist1;
  TH1D *mcHist1;
  TH1D *dataHist2;
  TH1D *mcHist2;

  TCanvas *c_data0 = (TCanvas*)data_file.Get("hist_trk_x_hit_Q_plane_0_prof");
  TCanvas *c_mc0   = (TCanvas*)mc_file.Get("hist_trk_x_hit_Q_plane_0_prof");
  dataHist0 = (TH1D*)c_data0->GetPrimitive("histtm_trk_x_hit_Q_plane_0");
  mcHist0 = (TH1D*)c_mc0->GetPrimitive("histtm_trk_x_hit_Q_plane_0");

  TCanvas *c_data1 = (TCanvas*)data_file.Get("hist_trk_x_hit_Q_plane_1_prof");
  TCanvas *c_mc1   = (TCanvas*)mc_file.Get("hist_trk_x_hit_Q_plane_1_prof");
  dataHist1 = (TH1D*)c_data1->GetPrimitive("histtm_trk_x_hit_Q_plane_1");
  mcHist1 = (TH1D*)c_mc1->GetPrimitive("histtm_trk_x_hit_Q_plane_1");

  TCanvas *c_data2 = (TCanvas*)data_file.Get("hist_trk_x_hit_Q_plane_2_prof");
  TCanvas *c_mc2   = (TCanvas*)mc_file.Get("hist_trk_x_hit_Q_plane_2_prof");
  dataHist2 = (TH1D*)c_data2->GetPrimitive("histtm_trk_x_hit_Q_plane_2");
  mcHist2 = (TH1D*)c_mc2->GetPrimitive("histtm_trk_x_hit_Q_plane_2");

  TCanvas* c = new TCanvas("c");
  mcHist0->Divide(dataHist0);
  mcHist0->SetLineColor(kBlue);
  mcHist0->SetLineWidth(2);
  mcHist1->Divide(dataHist1);
  mcHist1->SetLineColor(kRed);
  mcHist1->SetLineWidth(2);
  mcHist2->Divide(dataHist2);
  mcHist2->SetLineColor(kBlack);
  mcHist2->SetLineWidth(2);
  mcHist0->GetYaxis()->SetTitle("Hit Charge MC/Data Ratio");
  mcHist0->GetXaxis()->SetTitle("Track x Position [cm]");
  mcHist0->Draw();
  mcHist1->Draw("same");
  mcHist2->Draw("same");

  gStyle->SetOptStat(0);
  TLegend* leg = new TLegend(0.4,0.6,0.8,.8);  
  leg->AddEntry(mcHist0,"Plane 0", "l");
  leg->AddEntry(mcHist1,"Plane 1", "l");
  leg->AddEntry(mcHist2,"Plane 2", "l");
  leg->Draw("same");
		

}
