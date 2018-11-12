#include "fitData.cpp"

TTree *buildTree(unsigned int nEv, TF1 *fun) {
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

void toyMC2() {
  // gROOT->ProcessLine(".L fitData.cpp");
  gRandom = new TRandom3(0);

  float tau = 2.2e-6;
  auto fun = new TF1("exp", "[0] * exp(-x/[1])", 0, 1e-4);
  fun->SetParameters(1 / tau, tau);
  fun->SetNpx(10000);

  TGraph* tauAn = new TGraph();
  TGraph* tauAnRes = new TGraph();
  TGraph* sigmaAn = new TGraph();

  tauAn->SetTitle(";N;#tau (s)");
  sigmaAn->SetTitle(";N;#sigma (s)");

  std::vector<unsigned int> samples;
  unsigned int startN = 50, endN = 1000000, nN = 100;
  float dlogN = log10(endN / startN) / nN;

  for (float lN = log10(startN); lN <= log10(endN); lN += dlogN) {
    samples.push_back(TMath::Nint(pow(10, lN)));
  }

  for (unsigned int N : samples) {
    auto tree = buildTree(N, fun);

    std::map<std::string, measurement> results = fitData(tree);

    tauAn->SetPoint( tauAn->GetN(), N, results["analytical"].value );
    tauAnRes->SetPoint( tauAnRes->GetN(), N, std::fabs(results["analytical"].value - tau) );
    sigmaAn->SetPoint( sigmaAn->GetN(), N, results["analytical"].GetError() );

    delete tree;
  }

  auto trueSigma = new TF1("trueSigma", "[0] + [1]/sqrt(x)", 1, 1e7);
  trueSigma->SetLineColor(4);
  auto tS1 = new TF1("tS1", "[0] + [1]/sqrt(x)", 1, 1e7);
  tS1->SetLineColor(4);
  auto tS2 = new TF1("tS2", "[0] + [1]/sqrt(x)", 1, 1e7);
  tS2->SetLineColor(4);

  auto cc = new TCanvas("cc", "", 0, 0, 1600 + 4, 800 + 28);
  cc->Divide(2, 1);

  tauAn->SetMarkerStyle(20);
  tauAnRes->SetMarkerStyle(20);
  sigmaAn->SetMarkerStyle(20);
  tauAn->SetMarkerColor(2);
  tauAnRes->SetMarkerColor(2);
  sigmaAn->SetMarkerColor(2);

  cc->cd(1);
  tauAn->Draw("AP");
  tS1->SetParameter(0, tau);
  tS1->SetParameter(1, tau);
  tS1->Draw("same");
  tS2->SetParameter(0, tau);
  tS2->SetParameter(1, -tau);
  tS2->Draw("same");
  gPad->SetLogx();
  cc->cd(2);
  trueSigma->SetParameter(0, 0);
  trueSigma->SetParameter(1, tau);
  sigmaAn->Draw("AP");
  trueSigma->Draw("same");
  gPad->SetLogx();
  // cc->cd(3);
  // tauAnRes->Draw("AP");
  // gPad->SetLogx();
  // gPad->SetLogy();

  cc->Print("scaling.pdf");

}
