// /pnfs/uboone/data/uboone/reconstructed/prod_v08_00_00_18/data_extbnb_mcc9.1_v08_00_00_18/run1_reco2_C1/00/00/69/58/PhysicsRun-2016_7_24_11_54_40-0006958-00040_20160724T214243_ext_bnb_20160724T234647_merged_20181107T121731_optfilter_20181224T075125_reco1_postwcct_postdl_20181224T082906_20190723T182808_reco2.root


/// Psuedo code:
/*

  Find all the T0 from acpttrigtagger
  
  then grab the tracks that have an association through acpttrigtagger
  
  using these tracks get the associated caloritmetry data product

  then got through all the calo's trajectory points, mark it's x, then find the associated hit and store it's attributes 
  

 */


// Standard things to include
#include <iostream>
#include <stdlib.h>
#include <string>
#include <vector>

#include <fstream>
#include <iterator>
#include <algorithm>
#include <math.h> 
// These are the includes to use "Root" things 
#include "TInterpreter.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"

// These are the larsoft includes that let you
// have access to data-products and the event 
// details
#include "canvas/Utilities/InputTag.h"
#include "gallery/Event.h"
//#include "gallery/Event.h"
#include "gallery/ValidHandle.h"
#include "gallery/Handle.h"
#include "canvas/Persistency/Common/FindMany.h"
#include "canvas/Persistency/Common/FindOne.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/fwd.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Provenance/Timestamp.h"

//I'll need, calo, tracks, hits, anab::T0
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "canvas/Persistency/Provenance/EventAuxiliary.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/RawData/OpDetWaveform.h"

//This way you can be lazy
using namespace art;
using namespace std;

void Tracks_TreeMakerMC(){
  /*
  // create a vector of files we want to process
  std::vector<std::string> filenames;
  
  // read in a file list that we get from a sam-def but remember it 
  // is very long so if we want to run over it all it'll take a while
  // and we'll probably want to break it up on the grid
  ifstream myfile("file_lists/prod_crossingmuon_uboone_overlay_run1_reco2_reco2.list");
  copy(istream_iterator<string>(myfile),
       istream_iterator<string>(),
       back_inserter(filenames));
 
  //We'll just check the first 10k files for now (that's probably ~25k tracks)
  //filenames.erase(filenames.begin()+200,filenames.end());
  */
  
  //Here I just hard coding a file for testing, we can adjust this later 
  vector<string> filenames {"test.root"};
  
  

  // Here we will create all of our histograms 
  // I did this crazy inefficiently but I don't really care
  // This is currently only set up for single dimensional 
  // projections but extenting it to 2D will be straight forward
  TFile* out = new TFile("output.root","RECREATE");  
  TTree* fTree = new TTree("hit_v_track","HitPropertiesTrackProperties");
  int Run;
  fTree->Branch("Run",&Run);
  int Subrun;
  fTree->Branch("Subrun",&Subrun);
  int Event;
  fTree->Branch("Event",&Event);
  int Ntrks;
  fTree->Branch("Ntrks",&Ntrks);
  int Nacpt;
  fTree->Branch("Nacpt",&Nacpt);  
  int n_acpt;
  fTree->Branch("n_acpt",&n_acpt);  
  int Nhits;
  fTree->Branch("Nhits",&Nhits);  
  int n_hits;
  fTree->Branch("n_hits",&n_hits);  
  double hit_Q;
  double hit_A;
  double hit_sigma;
  double hit_time;
  int hit_plane;
  int hit_isMC;
  fTree->Branch("hit_Q",&hit_Q);  
  fTree->Branch("hit_A",&hit_A);  
  fTree->Branch("hit_sigma",&hit_sigma);  
  fTree->Branch("hit_time",&hit_time);  
  fTree->Branch("hit_plane",&hit_plane);  
  fTree->Branch("hit_isMC",&hit_isMC);  

  double trk_x;
  double trk_y;
  double trk_z;
  double trk_L;
  double trk_StartX;
  double trk_StartY;
  double trk_StartZ;
  double trk_EndX;
  double trk_EndY;
  double trk_EndZ;
  double trk_cosTheta;
  double trk_phi;
  double trk_ThetaXZ;
  double trk_ThetaYZ;
  double trk_dEdx;
  double trk_fracMC;

  fTree->Branch("trk_x",&trk_x);  
  fTree->Branch("trk_y",&trk_y);  
  fTree->Branch("trk_z",&trk_z);  
  fTree->Branch("trk_L",&trk_L);  

  fTree->Branch("trk_StartX",&trk_StartX);  
  fTree->Branch("trk_StartY",&trk_StartY);  
  fTree->Branch("trk_StartZ",&trk_StartZ);  

  fTree->Branch("trk_EndX",&trk_EndX);  
  fTree->Branch("trk_EndY",&trk_EndY);  
  fTree->Branch("trk_EndZ",&trk_EndZ);  

  fTree->Branch("trk_cosTheta",&trk_cosTheta);  
  fTree->Branch("trk_phi",&trk_phi);  
  fTree->Branch("trk_ThetaXZ",&trk_ThetaXZ);  
  fTree->Branch("trk_ThetaYZ",&trk_ThetaYZ);  
  fTree->Branch("trk_dEdx",&trk_dEdx);  
  fTree->Branch("trk_fracMC",&trk_fracMC);  

  double clus_pand_width;
  double clus_pand_start_tick;
  double clus_pand_end_tick;

  double clus_patrec_width;
  double clus_patrec_start_tick;
  double clus_patrec_end_tick;

  fTree->Branch("clus_pand_width",&clus_pand_width);  
  fTree->Branch("clus_pand_start_tick",&clus_pand_start_tick);  
  fTree->Branch("clus_pand_end_tick",&clus_pand_end_tick);  
  
  fTree->Branch("clus_patrec_width",&clus_patrec_width);  
  fTree->Branch("clus_patrec_start_tick",&clus_patrec_start_tick);  
  fTree->Branch("clus_patrec_end_tick",&clus_patrec_end_tick);  

  double MC_optIntegral;
  double Data_optIntegral;
  double Overlay_optIntegral;

  fTree->Branch("MC_optIntegral",&MC_optIntegral);  
  fTree->Branch("Data_optIntegral",&Data_optIntegral);  
  fTree->Branch("Overlay_optIntegral",&Overlay_optIntegral);  

  //First things fist we need to iterate through each event
  // gallery makes it easy to just hand a vector of files
  for (gallery::Event ev(filenames) ; !ev.atEnd(); ev.next()) {
    
    /// Prep our Branches
    Run = 0; //
    Subrun = 0; //
    Event = 0; //
    Ntrks = 0; //
    Nacpt = 0;//
    n_acpt = -1;//
    Nhits = 0;//
    n_hits = -1;//
    hit_Q = 0; //
    hit_A = 0; //
    hit_sigma = 0; //
    hit_time = 0; //
    hit_plane = 0; //
    trk_x = 0;//
    trk_y = 0;//
    trk_z = 0;//
    trk_cosTheta = 0;//
    trk_phi = 0;//
    trk_ThetaXZ = 0;//
    trk_ThetaYZ = 0;//
    trk_dEdx = 0;//

    trk_StartX = 0;//
    trk_StartY = 0;//
    trk_StartZ = 0;//
    trk_EndX = 0;//
    trk_EndY = 0;//
    trk_EndZ = 0; //
    
    
    trk_L = 0;//
    trk_fracMC = 0;//
    hit_isMC = 0;//

    clus_pand_width = 0;
    clus_pand_start_tick = 0;
    clus_pand_end_tick = 0;

    clus_patrec_width = 0;
    clus_patrec_start_tick = 0;
    clus_patrec_end_tick = 0;

    MC_optIntegral = 0;
    Data_optIntegral = 0;
    Overlay_optIntegral = 0;

    

    Run = ev.eventAuxiliary().run();
    Subrun = ev.eventAuxiliary().subRun();
    Event = ev.eventAuxiliary().event();
    
       
    //We will now found the events that have a ACPT in-time track
    auto const& t0s = *ev.getValidHandle<vector<anab::T0>>("acpttrigtagger");
    
    //Skipping those that don't
    if(t0s.size() == 0) continue;
    
    Nacpt = t0s.size();
    
    //This associates the T0 tag with the track
    auto const &t0_assoc_handle =
      ev.getValidHandle<art::Assns<anab::T0, recob::Track>>("acpttrigtagger");

    //Make a vector that will hold our tagged tracks
    std::vector<recob::Track> ACPT_tracks;
    
    // store the tagged tracks into that vector
    for(auto &ass : *t0_assoc_handle){
      art::Ptr<recob::Track> temp_trk = ass.second;
      ACPT_tracks.emplace_back(*temp_trk);
    }
    
    // Now we'll need to set things up to collect the calorimetry data
    // Start by saying which tracks we want:
    auto const & track_list_ptr = ev.getValidHandle<std::vector <recob::Track> >("pandora");
    auto const & track_list = (*track_list_ptr);
    
    Ntrks = track_list.size();

    //This let's us find which hits are assocaited to a give track trajectory point
    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(track_list_ptr, ev, "pandora"); 

    
    //Find all the hits that are matched 
    // to MCParticles
    //auto const& hit_handle = ev.getValidHandle<std::vector<recob::Hit>>("gaushit");
    //art::FindManyP<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hit_handle,ev,"gaushitTruthMatch");

    auto const &hit_mcpart_assn =
      *ev.getValidHandle<art::Assns<recob::Hit,simb::MCParticle,anab::BackTrackerHitMatchingData>>("gaushitTruthMatch");

    std::unordered_set<size_t> my_hit_key_mc;

    for(auto &pair : hit_mcpart_assn){
      my_hit_key_mc.insert(pair.first.key());
    }


    // Find the clusters that'll match our hits 
    //    art::FindManyP<recob::Cluster> patrec_clus_per_hit(hit_handle, ev, "pandoraPatRec:allOutcomes");   
    //art::FindManyP<recob::Cluster> pand_clus_per_hit(hit_handle, ev, "pandora");   

    //Let's loop through our tracks and find our calorimetry things
    for(auto &trk : ACPT_tracks){      
      n_acpt++;
      for(int itrk = 0; itrk < int(track_list.size()); itrk++){	
	if(trk.ID() == track_list.at(itrk).ID()){ 
	  
	  // Now we have a track which is matched to a ACPT crossing track
	  // now we want the hits for this track

	  // This is the vector of hits
	  auto vhit = fmthm.at(itrk);
	  Nhits = vhit.size();

	  // Check fraction that are truth:
	  int n_truthhits = 0;
	  for(int ht = 0; ht < Nhits; ht++){
	    //auto part_vec = particles_per_hit.at(vhit[ht].key());
	    //if(part_vec.size() > 0){
	    if(my_hit_key_mc.count(vhit[ht].key()) > 0){
	      n_truthhits++;
	    }
	  }

	  if(Nhits > 0){
	    trk_fracMC = double(n_truthhits)/double(Nhits);
	  }
	  else{
	    trk_fracMC = 0.;
	  }

	  if(trk_fracMC < 0.5) std::cout << "\t \t \t \t :::::::::: FAKE TRACK Run " << Run << " event : " << Event << std::endl; 
	  //	  if(trk_fracMC > 0.5) continue;

	  //Next we'll get the associated calorimetries things:
	  art::FindMany<anab::Calorimetry>  fmcal(track_list_ptr, ev, "pandoracaliSCE");//"pandoracali");
	  
 	  trk_L = trk.Length();
	  trk_StartX = trk.Start().X();
	  trk_StartY = trk.Start().Y();
	  trk_StartZ = trk.Start().Z();
	  trk_EndX = trk.End().X();
	  trk_EndY = trk.End().Y();
	  trk_EndZ = trk.End().Z();


	  // Grabbing the various optical waveforms
	  // Data Part  
	  //auto const & data_opt_list   = *ev.getValidHandle<std::vector <raw::OpDetWaveform> >("pmtreadout:OpdetBeamHighGain");
	  //	  auto const & mc_opt_list     = *ev.getValidHandle<std::vector <raw::OpDetWaveform> >("pmtreadoutnonoise:OpdetBeamHighGain");
	  // auto const & overlay_opt_list= *ev.getValidHandle<std::vector <raw::OpDetWaveform> >("mixer:OpdetBeamHighGain");

	  //Get Waveform Info
	  MC_optIntegral = 0;
	  Data_optIntegral = 0;
	  Overlay_optIntegral = 0;
	  // Iterates through the PMTs 
	  /*
	  for(int pmt = 0; pmt < data_opt_list.size(); pmt++){
	    if(data_opt_list[pmt].ChannelNumber() > 32) continue;
	    for(int tck = 0; tck < 1500; tck++){
	      Data_optIntegral   += data_opt_list[pmt][tck] - 2048;
	    }		
	  }
	  std::cout << "Data wf done " << std::endl;
	  for(int pmt = 0; pmt < mc_opt_list.size(); pmt++){
	    if(mc_opt_list[pmt].ChannelNumber() > 32) continue;
	    for(int tck = 0; tck < 1500; tck++){		    
	      MC_optIntegral     += mc_opt_list[pmt][tck] - 2048;
	    }
	  }
	  for(int pmt = 0; pmt < overlay_opt_list.size(); pmt++){
	    if(overlay_opt_list[pmt].ChannelNumber() > 32) continue;
	    for(int tck = 0; tck < 1500; tck++){		    
	      Overlay_optIntegral+= overlay_opt_list[pmt][tck] - 2048;
	    }
	  }
	  */
	  
	  // This is the vector of traj point info
	  // the Index() of this is the traj point of this track
	  auto vmeta = fmthm.data(itrk);
	  
	  // Now we can get the calorimetry points
	  std::vector<const anab::Calorimetry*> calos = fmcal.at(itrk);

	  // this will count which calo point we're on
	  int count = 0;

	  //iterate through the planes:
	  for(int pl = 0; pl < 3; pl++){
	    	 
	    //iterate through track meta points :	    
	    for(int vp = 0; vp < vmeta.size(); vp++){
	      
	      // store the track trajectory point index
	      // for the track meta point
	      int ind = vmeta[vp]->Index();

	      // check that the traj point is in the calorimetry point
	      // and belongs to the plane we are interested in 
	      if(track_list.at(itrk).HasValidPoint(ind) && vhit[vp]->WireID().Plane == pl){
		
		n_hits++;

       		//auto part_vec = particles_per_hit.at(vhit[vp].key());		
		//if(part_vec.size() > 0){
		if(my_hit_key_mc.count(vhit[vp].key()) > 0){
		  hit_isMC = 1;
		}

		// Get Cluster Info
		/*
		auto patrec_clus = patrec_clus_per_hit.at(vhit[vp].key());
		auto pand_clus = pand_clus_per_hit.at(vhit[vp].key());
				
		clus_pand_start_tick = pand_clus[0]->StartTick();
		clus_pand_end_tick = pand_clus[0]->EndTick();
		clus_pand_width = max(clus_pand_start_tick,clus_pand_end_tick) - min(clus_pand_start_tick,clus_pand_end_tick);		
		
		clus_patrec_start_tick = pand_clus[0]->StartTick();		  
		clus_patrec_end_tick = pand_clus[0]->EndTick();

		for(int pr = 0; pr < patrec_clus.size(); pr++){
		  
		  if(patrec_clus[pr]->StartTick() != clus_patrec_start_tick) 
		    clus_patrec_start_tick = patrec_clus[pr]->StartTick();

		  if(patrec_clus[pr]->EndTick() != clus_patrec_end_tick) 
		    clus_patrec_end_tick = patrec_clus[pr]->EndTick();

		}
	       	       
		clus_patrec_width = max(clus_patrec_start_tick,clus_patrec_end_tick) - min(clus_patrec_start_tick,clus_patrec_end_tick);
		*/

		// Grab the track traj point
		// WE DON'T CURRENTLY USE THIS
		// I kept it for testing purposes 
		auto trjp = track_list.at(itrk).TrajectoryPoint(ind);
		
		// Grab the calo point
		auto calp = calos[pl]->XYZ()[count];
		auto caldEdx = calos[pl]->dEdx()[count];
		
		// We need to calculate the angles 
		// of the calo points 
		double Phi = 0;
		double cosTheta = 0;
		double ThetaXZ = 0;
		double ThetaYZ = 0;
		
		if(count < vmeta.size()-1){
		  
		  auto angle = (calos[pl]->XYZ()[count]) - (calos[pl]->XYZ()[count+1]);  
		  
		  Phi = angle.Phi(); 
		  cosTheta = cos(angle.Theta());
		  ThetaXZ = atan2(angle.X(),angle.Z());
		  ThetaYZ = atan2(angle.Y(),angle.Z());		  

		}
		else{
		  auto angle = (calos[pl]->XYZ()[count-1]) - (calos[pl]->XYZ()[count]);  
		  
		  Phi = angle.Phi(); 
		  cosTheta = cos(angle.Theta());
		  ThetaXZ = atan2(angle.X(),angle.Z());
		  ThetaYZ = atan2(angle.Y(),angle.Z());
		}	

		// Grab the matched hit 
		auto hit = vhit[vp];


		hit_A = hit->PeakAmplitude();
		hit_Q = hit->Integral();
		hit_sigma = hit->RMS();
		hit_time = hit->PeakTime();
		hit_plane = pl;
		
		trk_x = calp.X(); 
		trk_y = calp.Y(); 
		trk_z = calp.Z(); 
		trk_phi = Phi;
		trk_cosTheta = cosTheta;
		trk_dEdx = caldEdx;
		trk_ThetaXZ = ThetaXZ;
		trk_ThetaYZ = ThetaYZ;
		
		fTree->Fill();
		
		// this tracks the correct calorimetry 
		// we are supposed to be anlayzing 
		count++;
	      }
	    }
	    count = 0; 
	  }
	}
	

      }
    }
  }

  out->cd();  
  fTree->Write();
  

}
  

