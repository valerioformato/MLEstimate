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
measurement fitMinuit(TTree *tree, std::string version);

void DivideByBinWidth(TH1 *hist) {
  for (int ibin = 0; ibin < hist->GetNbinsX(); ibin++) {
    float c = hist->GetBinContent(ibin + 1);
    float ec = hist->GetBinError(ibin + 1);
    float bw = hist->GetBinWidth(ibin + 1);
    hist->SetBinContent(ibin + 1, c / bw);
    hist->SetBinError(ibin + 1, ec / bw);
  }
};

std::map<std::string, measurement> fitData(TTree* tree) {

  std::map<std::string, measurement> map;

  float t;
  tree->SetBranchAddress("time", &t);
  Long64_t nEv = tree->GetEntries();
  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    tree->GetEntry(iEv);
  }

  map["analytical"] = fitAnalytical(tree);
  map["LL"] = fitScanLL(tree);
  map["Minuit"] = fitMinuit(tree, "Minuit");
  map["Minuit2"] = fitMinuit(tree, "Minuit2");

  return map;
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

  return tauML;
}

measurement fitScanLL(TTree *tree) {
  TCanvas *_canvas = (TCanvas *)gROOT->FindObject("ccL");

  float t;
  tree->SetBranchAddress("time", &t);
  Long64_t nEv = tree->GetEntries();
  std::vector<float> tvec(nEv);
  for (Long64_t iEv = 0; iEv < nEv; iEv++) {
    tree->GetEntry(iEv);
    tvec[iEv] = t;
  }

  measurement tauLL;
  tauLL.name = "tau";

  TGraph *logL = new TGraph();

  auto lfun = [](float _x, float _tau) { return 1. / _tau * exp(-_x / _tau); };
  double ll = 0;
  float tauStart = 2e-6, tauEnd = 2.8e-6, dtau = 1e-10;
  for (float tauThis = tauStart; tauThis <= tauEnd; tauThis += dtau) {
    ll = 0;
    for( auto tIt : tvec ){
      ll += log(lfun(tIt, tauThis));
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

  tauLL.variance = pow(0.5 * (t2 - t1), 2);

  cout << "LogLikelihood scan estimate" << endl;
  tauLL.Print();

  return tauLL;
}

measurement fitMinuit(TTree *tree, std::string method) {
  //create wrapper class for lilelihood
  class LL {
  public:
    LL(TTree *tree) : _tree(tree){ Init(); };

    void SetTree(TTree *tree) { _tree = tree; Init(); };

    void Init(){
      float t;
      _tree->SetBranchAddress("time", &t);
      Long64_t nEv = _tree->GetEntries();
      tvec.reserve(nEv);
      for (Long64_t iEv = 0; iEv < nEv; iEv++) {
        _tree->GetEntry(iEv);
        tvec.push_back(t);
      }
    }

    //this is particularly important
    double operator()(const double *_tau) {
      double ll = 0;
      for( auto tIt : tvec ){
        ll -= log(1. / *_tau * exp(-tIt / *_tau));
      }
      return ll;
    }

  private:
    TTree *_tree;
    std::vector<float> tvec;
  };

  auto lfun = new LL(tree);
  Long64_t nEv = tree->GetEntries();

  measurement tauMinuit;
  tauMinuit.name = "tau";

  auto minimizer = ROOT::Math::Factory::CreateMinimizer(method, "Migrad");
  minimizer->SetMaxFunctionCalls(1000000);
  minimizer->SetTolerance(0.001);
  minimizer->SetErrorDef(0.5);
  minimizer->SetPrintLevel(1);

  ROOT::Math::Functor fcn(*lfun, 1);
  double step[1] = {1e-9};
  double _tau[1] = {1e-6};

  minimizer->SetFunction(fcn);
  minimizer->SetVariable(0, "tau", _tau[0], step[0]);

  minimizer->Minimize();

  tauMinuit.value = minimizer->X()[0];
  tauMinuit.variance = minimizer->CovMatrix(0, 0);

  cout << method << " estimate" << endl;
  tauMinuit.Print();

  return tauMinuit;
}
