#include "fitData.cpp"

TTree *buildTree(unsigned int nEv, TF1* fun) {
  auto tree = new TTree("Triggers", "Triggers");

  float t = 0, dt = 0;
  tree->Branch("time", &t);

  cout << "Generating " << nEv << " events" << endl;

  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    t = fun->GetRandom();
    tree->Fill();
  }
  return tree;
}

void toyMC() {
  // gROOT->ProcessLine(".L fitData.cpp");
  gRandom = new TRandom3(0);

  float tau = 2.2e-6;
  auto fun = new TF1("exp", "[0] * exp(-x/[1])", 0, 1e-4);
  fun->SetParameters(1 / tau, tau);
  fun->SetNpx(10000);

  unsigned int nExp = 1000;
  unsigned int N = 475;

  TH1F *tauAn = new TH1F("tauAn", ";#tau (s);", 100, 2e-6, 2.6e-6);
  TH1F *tauLL = new TH1F("tauLL", ";#tau (s);", 100, 2e-6, 2.6e-6);
  TH1F *tauM = new TH1F("tauM", ";#tau (s);", 100, 2e-6, 2.6e-6);
  TH1F *tauM2 = new TH1F("tauM2", ";#tau (s);", 100, 2e-6, 2.6e-6);

  TH1F *sigmaAn = new TH1F("sigmaAn", ";#sigma (s);", 100, 0.5e-7, 2.e-7);
  TH1F *sigmaLL = new TH1F("sigmaLL", ";#sigma (s);", 100, 0.5e-7, 2.e-7);
  TH1F *sigmaM = new TH1F("sigmaM", ";#sigma (s);", 100, 0.5e-7, 2.e-7);
  TH1F *sigmaM2 = new TH1F("sigmaM2", ";#sigma (s);", 100, 0.5e-7, 2.e-7);


  for (int iExp = 0; iExp < nExp; iExp++) {
    auto tree = buildTree(N, fun);

    std::map<std::string, measurement> results = fitData(tree);

    tauAn->Fill( results["analytical"].value );
    tauLL->Fill( results["LL"].value );
    tauM->Fill( results["Minuit"].value );
    tauM2->Fill( results["Minuit2"].value );

    sigmaAn->Fill( results["analytical"].GetError() );
    sigmaLL->Fill( results["LL"].GetError() );
    sigmaM->Fill( results["Minuit"].GetError() );
    sigmaM2->Fill( results["Minuit2"].GetError() );

    delete tree;
  }

  auto cc = new TCanvas("cc", "", 0, 0, 1600+4, 800+28);
  cc->Divide(4, 2);

  cc->cd(1);
  tauAn->Draw();
  cc->cd(2);
  tauLL->Draw();
  cc->cd(3);
  tauM->Draw();
  cc->cd(4);
  tauM2->Draw();
  cc->cd(5);
  sigmaAn->Draw();
  cc->cd(6);
  sigmaLL->Draw();
  cc->cd(7);
  sigmaM->Draw();
  cc->cd(8);
  sigmaM2->Draw();

}
