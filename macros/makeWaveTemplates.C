{
  //320mu 3-4 X0 - t_max[1] used as time reference
  TChain *c= new TChain("H4treeReco");
  c->Add("/tmp/meridian/RECO_3350.root");
  c->Add("/tmp/meridian/RECO_3351.root");
  c->Add("/tmp/meridian/RECO_3352.root");
  c->Add("/tmp/meridian/RECO_3353.root");
  c->Add("/tmp/meridian/RECO_3356.root");
  TProfile *h0=new TProfile("MCPTrigger_waveProfile","MCPTrigger_waveProfile",280,-4,3);
  TProfile *h1=new TProfile("SiPad1_waveProfile","SiPad1_waveProfile",280,-9,-2);
  TProfile *h2=new TProfile("SiPad2_waveProfile","SiPad2_waveProfile",280,-9,-2);
  c->Draw("(wave_aroundmax[1]/wave_max[1]):(time_aroundmax[1]*1e9-t_max[1])>>MCPTrigger_waveProfile","wave_max[1]>300","PROF");
  c->Draw("(wave_aroundmax[2]/wave_max[2]):(time_aroundmax[2]*1e9-t_max[1])>>SiPad1_waveProfile","wave_max[2]>50 && abs(t_max[2]-45)<8 && wave_max[1]>100 && abs(t_max_frac50[2]-t_max_frac50[1]+5.5)<1","PROF");
  c->Draw("(wave_aroundmax[3]/wave_max[3]):(time_aroundmax[3]*1e9-t_max[1])>>SiPad2_waveProfile","wave_max[3]>50 && abs(t_max[3]-45)<8 && wave_max[1]>100 && abs(t_max_frac50[3]-t_max_frac50[1]+5.1)<1","PROF");
  TFile *f=new TFile("waveTemplates.root","RECREATE");
  h0->Write("MCPTrigger_waveProfile");
  h1->Write("SiPad1_waveProfile");
  h2->Write("SiPad2_waveProfile");
}
