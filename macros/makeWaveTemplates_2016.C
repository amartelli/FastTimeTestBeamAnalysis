{
  //120mu Ele E 100 - t_max[1] used as time reference
  TChain* c120 = new TChain("H4treeReco");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3729.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3733.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3737.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3751.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3752.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3753.root");
  //E150
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3719.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3725.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3713.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3714.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3717.root");
  TProfile *h1=new TProfile("SiPad1_waveProfile","SiPad1_waveProfile",400, -10, 10);
  TProfile *h2=new TProfile("SiPad2_waveProfile","SiPad2_waveProfile",400, -10, 10);
  TProfile *h3=new TProfile("SiPad3_waveProfile","SiPad3_waveProfile",400, -10, 10);
  TProfile *h4=new TProfile("SiPad4_waveProfile","SiPad4_waveProfile",400, -10, 10);
  TProfile *h5=new TProfile("SiPad5_waveProfile","SiPad5_waveProfile",400, -10, 10);
  TProfile *h6=new TProfile("SiPad6_waveProfile","SiPad6_waveProfile",400, -10, 10);
  c120->Draw("(wave_aroundmax[4]/wave_max[4]):(time_aroundmax[4]*1e9-t_max[4])>>SiPad1_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4])","PROF");
  c120->Draw("(wave_aroundmax[5]/wave_max[5]):(time_aroundmax[5]*1e9-t_max[5])>>SiPad2_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[5]>(pedestal[5]+5*pedestalRMS[5])","PROF");
  c120->Draw("(wave_aroundmax[6]/wave_max[6]):(time_aroundmax[6]*1e9-t_max[6])>>SiPad3_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[6]>(pedestal[6]+5*pedestalRMS[6])","PROF");
  c120->Draw("(wave_aroundmax[3]/wave_max[3]):(time_aroundmax[3]*1e9-t_max[3])>>SiPad4_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[3]>(pedestal[3]+5*pedestalRMS[3])","PROF");
  c120->Draw("(wave_aroundmax[2]/wave_max[2]):(time_aroundmax[2]*1e9-t_max[2])>>SiPad5_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[2]>(pedestal[2]+5*pedestalRMS[2])","PROF");
  c120->Draw("(wave_aroundmax[1]/wave_max[1]):(time_aroundmax[1]*1e9-t_max[1])>>SiPad6_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[1]>(pedestal[1]+5*pedestalRMS[1])","PROF");

  TFile *f=new TFile("waveTemplates2016_120um.root","RECREATE");
  h1->Write("SiPad1_waveProfile");
  h2->Write("SiPad2_waveProfile");
  h3->Write("SiPad3_waveProfile");
  h4->Write("SiPad4_waveProfile");
  h5->Write("SiPad5_waveProfile");
  h6->Write("SiPad6_waveProfile");
  f->Close();
  c120->Delete();
  h1->Delete();
  h2->Delete();
  h3->Delete();
  h4->Delete();
  h5->Delete();
  h6->Delete();

  //200um
  c120 = new TChain("H4treeReco");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3767.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3768.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3756.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3766.root");
  //E150
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3769.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3771.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3772.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3773.root");
  h1=new TProfile("SiPad1_waveProfile","SiPad1_waveProfile",400, -10, 10);
  h2=new TProfile("SiPad2_waveProfile","SiPad2_waveProfile",400, -10, 10);
  h3=new TProfile("SiPad3_waveProfile","SiPad3_waveProfile",400, -10, 10);
  h4=new TProfile("SiPad4_waveProfile","SiPad4_waveProfile",400, -10, 10);
  h5=new TProfile("SiPad5_waveProfile","SiPad5_waveProfile",400, -10, 10);
  h6=new TProfile("SiPad6_waveProfile","SiPad6_waveProfile",400, -10, 10);
  c120->Draw("(wave_aroundmax[4]/wave_max[4]):(time_aroundmax[4]*1e9-t_max[4])>>SiPad1_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4])","PROF");
  c120->Draw("(wave_aroundmax[5]/wave_max[5]):(time_aroundmax[5]*1e9-t_max[5])>>SiPad2_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[5]>(pedestal[5]+5*pedestalRMS[5])","PROF");
  c120->Draw("(wave_aroundmax[6]/wave_max[6]):(time_aroundmax[6]*1e9-t_max[6])>>SiPad3_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[6]>(pedestal[6]+5*pedestalRMS[6])","PROF");
  c120->Draw("(wave_aroundmax[3]/wave_max[3]):(time_aroundmax[3]*1e9-t_max[3])>>SiPad4_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[3]>(pedestal[3]+5*pedestalRMS[3])","PROF");
  c120->Draw("(wave_aroundmax[2]/wave_max[2]):(time_aroundmax[2]*1e9-t_max[2])>>SiPad5_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[2]>(pedestal[2]+5*pedestalRMS[2])","PROF");
  c120->Draw("(wave_aroundmax[1]/wave_max[1]):(time_aroundmax[1]*1e9-t_max[1])>>SiPad6_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[1]>(pedestal[1]+5*pedestalRMS[1])","PROF");
  f = new TFile("waveTemplates2016_200um.root","RECREATE");
  h1->Write("SiPad1_waveProfile");
  h2->Write("SiPad2_waveProfile");
  h3->Write("SiPad3_waveProfile");
  h4->Write("SiPad4_waveProfile");
  h5->Write("SiPad5_waveProfile");
  h6->Write("SiPad6_waveProfile");
  f->Close();
  c120->Delete();
  h1->Delete();
  h2->Delete();
  h3->Delete();
  h4->Delete();
  h5->Delete();
  h6->Delete();



  //300um  E 150
  c120 = new TChain("H4treeReco");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3777.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3778.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3774.root");
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3775.root");
  //E100
  c120->Add("root://eoscms//store/group/upgrade/HGCAL/TimingTB_H2_Apr2016/recoTest/RECO_3779.root");
  h1=new TProfile("SiPad1_waveProfile","SiPad1_waveProfile",400, -10, 10);
  h2=new TProfile("SiPad2_waveProfile","SiPad2_waveProfile",400, -10, 10);
  h3=new TProfile("SiPad3_waveProfile","SiPad3_waveProfile",400, -10, 10);
  h4=new TProfile("SiPad4_waveProfile","SiPad4_waveProfile",400, -10, 10);
  h5=new TProfile("SiPad5_waveProfile","SiPad5_waveProfile",400, -10, 10);
  h6=new TProfile("SiPad6_waveProfile","SiPad6_waveProfile",400, -10, 10);
  c120->Draw("(wave_aroundmax[4]/wave_max[4]):(time_aroundmax[4]*1e9-t_max[4])>>SiPad1_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4])","PROF");
  c120->Draw("(wave_aroundmax[5]/wave_max[5]):(time_aroundmax[5]*1e9-t_max[5])>>SiPad2_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[5]>(pedestal[5]+5*pedestalRMS[5])","PROF");
  c120->Draw("(wave_aroundmax[6]/wave_max[6]):(time_aroundmax[6]*1e9-t_max[6])>>SiPad3_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[6]>(pedestal[6]+5*pedestalRMS[6])","PROF");
  c120->Draw("(wave_aroundmax[3]/wave_max[3]):(time_aroundmax[3]*1e9-t_max[3])>>SiPad4_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[3]>(pedestal[3]+5*pedestalRMS[3])","PROF");
  c120->Draw("(wave_aroundmax[2]/wave_max[2]):(time_aroundmax[2]*1e9-t_max[2])>>SiPad5_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[2]>(pedestal[2]+5*pedestalRMS[2])","PROF");
  c120->Draw("(wave_aroundmax[1]/wave_max[1]):(time_aroundmax[1]*1e9-t_max[1])>>SiPad6_waveProfile","wave_max[4]>(pedestal[4]+5*pedestalRMS[4]) && wave_max[1]>(pedestal[1]+5*pedestalRMS[1])","PROF");
  f = new TFile("waveTemplates2016_300um.root","RECREATE");
  h1->Write("SiPad1_waveProfile");
  h2->Write("SiPad2_waveProfile");
  h3->Write("SiPad3_waveProfile");
  h4->Write("SiPad4_waveProfile");
  h5->Write("SiPad5_waveProfile");
  h6->Write("SiPad6_waveProfile");
  f->Close();
  c120->Delete();
  h1->Delete();
  h2->Delete();
  h3->Delete();
  h4->Delete();
  h5->Delete();
  h6->Delete();

  std::cout << " CIAO " << std::endl;


}
