float tau = 2.2e-6;

struct measurement {
  std::string name;
  float value = 0;
  float variance = 0;
  float GetError() { return sqrt(variance); };
  void Print() {
    cout << " - " << name << " = ";
    cout << value << " +- " << GetError();
    cout << "   (variance = " << variance << ")" << endl;
  }
};

measurement fitAnalytical(TTree *tree);
measurement fitScanLL(TTree *tree);
measurement fitMinuit(TTree *tree);
measurement fitMinuit2(TTree *tree);

void DivideByBinWidth(TH1 *hist) {
  for (int ibin = 0; ibin < hist->GetNbinsX(); ibin++) {
    float c = hist->GetBinContent(ibin + 1);
    float ec = hist->GetBinError(ibin + 1);
    float bw = hist->GetBinWidth(ibin + 1);
    hist->SetBinContent(ibin + 1, c / bw);
    hist->SetBinError(ibin + 1, ec / bw);
  }
};

auto infile = std::unique_ptr<TFile>(TFile::Open("data/data.root"));

void fitData() {

  auto hist = new TH1F("hist", ";t (s);Counts", 1000, 0, 1e-4);
  auto tree = (TTree *)infile->Get("Triggers");
  float t;
  tree->SetBranchAddress("time", &t);
  Long64_t nEv = tree->GetEntries();
  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    tree->GetEntry(iEv);
    hist->Fill(t);
  }

  DivideByBinWidth(hist);
  hist->SetStats(0);
  hist->SetMarkerStyle(20);
  hist->SetLineColor(kBlack);
  hist->GetXaxis()->SetRangeUser(0, 1.1e-5);

  auto fun = new TF1("fun", "[0] * exp(-x/[1])", 0, 1e-4);
  fun->SetParameters(nEv / tau, tau);
  fun->SetNpx(10000);

  auto ccA = new TCanvas("ccA", "Analytical solution");
  ccA->Draw();
  hist->Draw("E");
  fun->Draw("SAME");
  auto tauAn = fitAnalytical(tree);
  ccA->Print("plots/plot2_An.pdf");

  auto ccL = new TCanvas("ccL", "LogLikelihood scan");
  ccL->Draw();
  hist->Draw("E");
  fun->Draw("SAME");
  auto tauLL = fitScanLL(tree);
  ccA->Print("plots/plot2_LL.pdf");
}

measurement fitAnalytical(TTree *tree) {
  float t;
  tree->SetBranchAddress("time", &t);

  measurement tauML;
  tauML.name = "tau";

  Long64_t nEv = tree->GetEntries();
  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    tree->GetEntry(iEv);

    tauML.value += t;
  }
  tauML.value /= nEv;
  tauML.variance = (tauML.value * tauML.value) / nEv;

  cout << "Analytical estimate" << endl;
  tauML.Print();

  auto funA = new TF1("funA", "[0] * exp(-x/[1])", 0, 1e-4);
  funA->SetParameters(nEv / tauML.value, tauML.value);
  funA->SetLineColor(4);
  funA->SetNpx(10000);
  funA->Draw("SAME");

  return tauML;
}

measurement fitScanLL(TTree *tree) {
  TCanvas *_canvas = (TCanvas *)gROOT->FindObject("ccL");

  float t;
  tree->SetBranchAddress("time", &t);

  measurement tauLL;
  tauLL.name = "tau";

  Long64_t nEv = tree->GetEntries();

  TGraph *logL = new TGraph();

  auto lfun = [](float _x, float _tau) { return 1. / _tau * exp(-_x / _tau); };
  double ll = 0;
  float tauStart = 2e-6, tauEnd = 2.8e-6, dtau = 1e-10;
  for (float tauThis = tauStart; tauThis <= tauEnd; tauThis += dtau) {
    ll = 0;
    for (Long64_t iEv = 0; iEv < nEv; iEv++) {
      tree->GetEntry(iEv);

      ll += log(lfun(t, tauThis));
    }
    logL->SetPoint(logL->GetN(), tauThis, ll);
  }

  ll = logL->GetY()[0];
  tauLL.value = logL->GetX()[0];
  for (int ip = 0; ip < logL->GetN(); ip++) {
    if (logL->GetY()[ip] > ll) {
      ll = logL->GetY()[ip];
      tauLL.value = logL->GetX()[ip];
    }
  }
  float t1 = 0, t2 = 0;
  for (int ip = 0; ip < logL->GetN(); ip++) {
    if (logL->GetY()[ip] > ll - 0.5 && !(t1 > 0)) {
      t1 = logL->GetX()[ip];
    }
    if (logL->GetY()[ip] < ll - 0.5 && (t1 > 0) && !(t2 > 0)) {
      t2 = logL->GetX()[ip];
    }
  }
  // cout << t1 << " " << t2 << endl;

  tauLL.variance = pow(0.5 * (t2 - t1), 2);

  cout << "LogLikelihood scan estimate" << endl;
  tauLL.Print();

  _canvas->cd(1);
  auto funLL = new TF1("funLL", "[0] * exp(-x/[1])", 0, 1e-4);
  funLL->SetParameters(nEv / tauLL.value, tauLL.value);
  funLL->SetLineColor(4);
  funLL->SetNpx(10000);
  funLL->Draw("SAME");

  auto _cc = new TCanvas();
  float miny = std::min(logL->Eval(2.1e-6), logL->Eval(2.65e-6));
  logL->GetXaxis()->SetRangeUser(2.1e-6, 2.65e-6);
  logL->GetYaxis()->SetRangeUser(miny, ll + 1);
  logL->SetTitle(";#tau (s);logL");
  logL->SetLineWidth(2);
  logL->Draw("AL");

  auto l1 = new TLine(tauLL.value, miny, tauLL.value, ll);
  l1->SetLineColor(2);
  l1->Draw("SAME");

  auto l2 = new TLine(t1, ll-0.5, t2, ll-0.5);
  l2->SetLineColor(2);
  l2->SetLineStyle(2);
  l2->Draw("SAME");

  _cc->Print("plots/plot2_LLscan.pdf");

  return tauLL;
}
