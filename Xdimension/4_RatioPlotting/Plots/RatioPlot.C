void RatioPlot()
{

  TFile data_file("../out_hists_data.root","read");
  TFile mc_file("../out_hists_mc.root","read");
  TH1D *dataHist;
  TH1D *mcHist;
  TCanvas *c_data = (TCanvas*)data_file.Get("hist_trk_x_hit_sigma_plane_2_prof");
  TCanvas *c_mc   = (TCanvas*)mc_file.Get("hist_trk_x_hit_sigma_plane_2_prof");
  dataHist = (TH1D*)c_data->GetPrimitive("histtm_trk_x_hit_sigma_plane_2");
  mcHist = (TH1D*)c_mc->GetPrimitive("histtm_trk_x_hit_sigma_plane_2");

  TCanvas* c = new TCanvas("c");
  mcHist->Divide(dataHist);
  mcHist->GetYaxis()->SetTitle("MC/Data");
  mcHist->GetXaxis()->SetTitle("Track x Position [cm]");
  mcHist->Draw();

  
}
