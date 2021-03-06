{
//=========Macro generated from canvas: c1/c1
//=========  (Wed Jun  4 17:22:36 2014) by ROOT version5.34/14
   TCanvas *c1 = new TCanvas("c1", "c1",66,52,700,500);
   c1->Range(-1037.5,0.540949,4837.5,6.270892);
   c1->SetFillColor(0);
   c1->SetBorderMode(0);
   c1->SetBorderSize(2);
   c1->SetLogy();
   c1->SetFrameBorderMode(0);
   c1->SetFrameBorderMode(0);
   
   TH1F *hAll = new TH1F("hAll","ScintE",100,-450,4250);
   hAll->SetBinContent(9,267);
   hAll->SetBinContent(10,263238);
   hAll->SetBinContent(11,38371);
   hAll->SetBinContent(12,21075);
   hAll->SetBinContent(13,14478);
   hAll->SetBinContent(14,10507);
   hAll->SetBinContent(15,8080);
   hAll->SetBinContent(16,6383);
   hAll->SetBinContent(17,5185);
   hAll->SetBinContent(18,4348);
   hAll->SetBinContent(19,3574);
   hAll->SetBinContent(20,2996);
   hAll->SetBinContent(21,2406);
   hAll->SetBinContent(22,2161);
   hAll->SetBinContent(23,1864);
   hAll->SetBinContent(24,1607);
   hAll->SetBinContent(25,1333);
   hAll->SetBinContent(26,1181);
   hAll->SetBinContent(27,1011);
   hAll->SetBinContent(28,950);
   hAll->SetBinContent(29,782);
   hAll->SetBinContent(30,765);
   hAll->SetBinContent(31,629);
   hAll->SetBinContent(32,589);
   hAll->SetBinContent(33,516);
   hAll->SetBinContent(34,512);
   hAll->SetBinContent(35,447);
   hAll->SetBinContent(36,369);
   hAll->SetBinContent(37,359);
   hAll->SetBinContent(38,328);
   hAll->SetBinContent(39,344);
   hAll->SetBinContent(40,334);
   hAll->SetBinContent(41,314);
   hAll->SetBinContent(42,279);
   hAll->SetBinContent(43,318);
   hAll->SetBinContent(44,293);
   hAll->SetBinContent(45,255);
   hAll->SetBinContent(46,254);
   hAll->SetBinContent(47,233);
   hAll->SetBinContent(48,215);
   hAll->SetBinContent(49,187);
   hAll->SetBinContent(50,171);
   hAll->SetBinContent(51,184);
   hAll->SetBinContent(52,164);
   hAll->SetBinContent(53,159);
   hAll->SetBinContent(54,145);
   hAll->SetBinContent(55,136);
   hAll->SetBinContent(56,154);
   hAll->SetBinContent(57,146);
   hAll->SetBinContent(58,130);
   hAll->SetBinContent(59,147);
   hAll->SetBinContent(60,140);
   hAll->SetBinContent(61,134);
   hAll->SetBinContent(62,158);
   hAll->SetBinContent(63,156);
   hAll->SetBinContent(64,160);
   hAll->SetBinContent(65,169);
   hAll->SetBinContent(66,195);
   hAll->SetBinContent(67,224);
   hAll->SetBinContent(68,267);
   hAll->SetBinContent(69,351);
   hAll->SetBinContent(70,491);
   hAll->SetBinContent(71,575);
   hAll->SetBinContent(72,696);
   hAll->SetBinContent(73,860);
   hAll->SetBinContent(74,986);
   hAll->SetBinContent(75,1108);
   hAll->SetBinContent(76,1159);
   hAll->SetBinContent(77,1262);
   hAll->SetBinContent(78,1168);
   hAll->SetBinContent(79,1149);
   hAll->SetBinContent(80,1011);
   hAll->SetBinContent(81,887);
   hAll->SetBinContent(82,786);
   hAll->SetBinContent(83,611);
   hAll->SetBinContent(84,579);
   hAll->SetBinContent(85,443);
   hAll->SetBinContent(86,317);
   hAll->SetBinContent(87,258);
   hAll->SetBinContent(88,210);
   hAll->SetBinContent(89,134);
   hAll->SetBinContent(90,104);
   hAll->SetBinContent(91,67);
   hAll->SetBinContent(92,26);
   hAll->SetBinContent(93,3221);
   hAll->SetEntries(420935);
   
   TPaveStats *ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   TText *text = ptstats->AddText("hAll");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 420935 ");
   text = ptstats->AddText("Mean  =  247.2");
   text = ptstats->AddText("RMS   =  721.9");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   hAll->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(hAll);
   hAll->GetXaxis()->SetRange(1,100);
   hAll->GetXaxis()->SetLabelFont(42);
   hAll->GetXaxis()->SetLabelSize(0.035);
   hAll->GetXaxis()->SetTitleSize(0.035);
   hAll->GetXaxis()->SetTitleFont(42);
   hAll->GetYaxis()->SetLabelFont(42);
   hAll->GetYaxis()->SetLabelSize(0.035);
   hAll->GetYaxis()->SetTitleSize(0.035);
   hAll->GetYaxis()->SetTitleFont(42);
   hAll->GetZaxis()->SetLabelFont(42);
   hAll->GetZaxis()->SetLabelSize(0.035);
   hAll->GetZaxis()->SetTitleSize(0.035);
   hAll->GetZaxis()->SetTitleFont(42);
   hAll->Draw("");
   
   TPaveText *pt = new TPaveText(0.4368391,0.94,0.5631609,0.995,"blNDC");
   pt->SetName("title");
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->SetFillStyle(0);
   pt->SetTextFont(42);
   text = pt->AddText("ScintE");
   pt->Draw();
   
   TH1F *h0 = new TH1F("h0","ScintE {PID==0}",100,-450,4250);
   h0->SetBinContent(9,95);
   h0->SetBinContent(10,179211);
   h0->SetBinContent(11,32196);
   h0->SetBinContent(12,15625);
   h0->SetBinContent(13,9268);
   h0->SetBinContent(14,5763);
   h0->SetBinContent(15,3977);
   h0->SetBinContent(16,2960);
   h0->SetBinContent(17,2300);
   h0->SetBinContent(18,1980);
   h0->SetBinContent(19,1643);
   h0->SetBinContent(20,1391);
   h0->SetBinContent(21,1170);
   h0->SetBinContent(22,1047);
   h0->SetBinContent(23,926);
   h0->SetBinContent(24,761);
   h0->SetBinContent(25,673);
   h0->SetBinContent(26,609);
   h0->SetBinContent(27,512);
   h0->SetBinContent(28,497);
   h0->SetBinContent(29,377);
   h0->SetBinContent(30,384);
   h0->SetBinContent(31,327);
   h0->SetBinContent(32,300);
   h0->SetBinContent(33,278);
   h0->SetBinContent(34,285);
   h0->SetBinContent(35,253);
   h0->SetBinContent(36,213);
   h0->SetBinContent(37,214);
   h0->SetBinContent(38,162);
   h0->SetBinContent(39,215);
   h0->SetBinContent(40,215);
   h0->SetBinContent(41,217);
   h0->SetBinContent(42,132);
   h0->SetBinContent(43,59);
   h0->SetBinContent(44,69);
   h0->SetBinContent(45,54);
   h0->SetBinContent(46,61);
   h0->SetBinContent(47,55);
   h0->SetBinContent(48,56);
   h0->SetBinContent(49,46);
   h0->SetBinContent(50,40);
   h0->SetBinContent(51,41);
   h0->SetBinContent(52,41);
   h0->SetBinContent(53,53);
   h0->SetBinContent(54,36);
   h0->SetBinContent(55,41);
   h0->SetBinContent(56,37);
   h0->SetBinContent(57,31);
   h0->SetBinContent(58,33);
   h0->SetBinContent(59,42);
   h0->SetBinContent(60,35);
   h0->SetBinContent(61,37);
   h0->SetBinContent(62,34);
   h0->SetBinContent(63,34);
   h0->SetBinContent(64,35);
   h0->SetBinContent(65,25);
   h0->SetBinContent(66,33);
   h0->SetBinContent(67,23);
   h0->SetBinContent(68,25);
   h0->SetBinContent(69,26);
   h0->SetBinContent(70,24);
   h0->SetBinContent(71,18);
   h0->SetBinContent(72,20);
   h0->SetBinContent(73,22);
   h0->SetBinContent(74,22);
   h0->SetBinContent(75,24);
   h0->SetBinContent(76,22);
   h0->SetBinContent(77,19);
   h0->SetBinContent(78,19);
   h0->SetBinContent(79,16);
   h0->SetBinContent(80,15);
   h0->SetBinContent(81,13);
   h0->SetBinContent(82,23);
   h0->SetBinContent(83,17);
   h0->SetBinContent(84,16);
   h0->SetBinContent(85,14);
   h0->SetBinContent(86,14);
   h0->SetBinContent(87,15);
   h0->SetBinContent(88,17);
   h0->SetBinContent(89,11);
   h0->SetBinContent(90,10);
   h0->SetBinContent(91,9);
   h0->SetBinContent(92,1);
   h0->SetBinContent(93,583);
   h0->SetEntries(268242);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   text = ptstats->AddText("h0");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 268242 ");
   text = ptstats->AddText("Mean  =  87.92");
   text = ptstats->AddText("RMS   =  306.8");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   h0->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(h0);
   h0->SetLineColor(9);
   h0->GetXaxis()->SetRange(1,100);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetLabelSize(0.035);
   h0->GetXaxis()->SetTitleSize(0.035);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetLabelSize(0.035);
   h0->GetYaxis()->SetTitleSize(0.035);
   h0->GetYaxis()->SetTitleFont(42);
   h0->GetZaxis()->SetLabelFont(42);
   h0->GetZaxis()->SetLabelSize(0.035);
   h0->GetZaxis()->SetTitleSize(0.035);
   h0->GetZaxis()->SetTitleFont(42);
   h0->Draw("Same");
   
   TH1F *h0 = new TH1F("h0","ScintE {PID==0}",100,-450,4250);
   h0->SetBinContent(9,95);
   h0->SetBinContent(10,179211);
   h0->SetBinContent(11,32196);
   h0->SetBinContent(12,15625);
   h0->SetBinContent(13,9268);
   h0->SetBinContent(14,5763);
   h0->SetBinContent(15,3977);
   h0->SetBinContent(16,2960);
   h0->SetBinContent(17,2300);
   h0->SetBinContent(18,1980);
   h0->SetBinContent(19,1643);
   h0->SetBinContent(20,1391);
   h0->SetBinContent(21,1170);
   h0->SetBinContent(22,1047);
   h0->SetBinContent(23,926);
   h0->SetBinContent(24,761);
   h0->SetBinContent(25,673);
   h0->SetBinContent(26,609);
   h0->SetBinContent(27,512);
   h0->SetBinContent(28,497);
   h0->SetBinContent(29,377);
   h0->SetBinContent(30,384);
   h0->SetBinContent(31,327);
   h0->SetBinContent(32,300);
   h0->SetBinContent(33,278);
   h0->SetBinContent(34,285);
   h0->SetBinContent(35,253);
   h0->SetBinContent(36,213);
   h0->SetBinContent(37,214);
   h0->SetBinContent(38,162);
   h0->SetBinContent(39,215);
   h0->SetBinContent(40,215);
   h0->SetBinContent(41,217);
   h0->SetBinContent(42,132);
   h0->SetBinContent(43,59);
   h0->SetBinContent(44,69);
   h0->SetBinContent(45,54);
   h0->SetBinContent(46,61);
   h0->SetBinContent(47,55);
   h0->SetBinContent(48,56);
   h0->SetBinContent(49,46);
   h0->SetBinContent(50,40);
   h0->SetBinContent(51,41);
   h0->SetBinContent(52,41);
   h0->SetBinContent(53,53);
   h0->SetBinContent(54,36);
   h0->SetBinContent(55,41);
   h0->SetBinContent(56,37);
   h0->SetBinContent(57,31);
   h0->SetBinContent(58,33);
   h0->SetBinContent(59,42);
   h0->SetBinContent(60,35);
   h0->SetBinContent(61,37);
   h0->SetBinContent(62,34);
   h0->SetBinContent(63,34);
   h0->SetBinContent(64,35);
   h0->SetBinContent(65,25);
   h0->SetBinContent(66,33);
   h0->SetBinContent(67,23);
   h0->SetBinContent(68,25);
   h0->SetBinContent(69,26);
   h0->SetBinContent(70,24);
   h0->SetBinContent(71,18);
   h0->SetBinContent(72,20);
   h0->SetBinContent(73,22);
   h0->SetBinContent(74,22);
   h0->SetBinContent(75,24);
   h0->SetBinContent(76,22);
   h0->SetBinContent(77,19);
   h0->SetBinContent(78,19);
   h0->SetBinContent(79,16);
   h0->SetBinContent(80,15);
   h0->SetBinContent(81,13);
   h0->SetBinContent(82,23);
   h0->SetBinContent(83,17);
   h0->SetBinContent(84,16);
   h0->SetBinContent(85,14);
   h0->SetBinContent(86,14);
   h0->SetBinContent(87,15);
   h0->SetBinContent(88,17);
   h0->SetBinContent(89,11);
   h0->SetBinContent(90,10);
   h0->SetBinContent(91,9);
   h0->SetBinContent(92,1);
   h0->SetBinContent(93,583);
   h0->SetEntries(268242);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   text = ptstats->AddText("h0");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 268242 ");
   text = ptstats->AddText("Mean  =  87.92");
   text = ptstats->AddText("RMS   =  306.8");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   h0->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(h0);
   h0->SetLineColor(9);
   h0->GetXaxis()->SetRange(1,100);
   h0->GetXaxis()->SetLabelFont(42);
   h0->GetXaxis()->SetLabelSize(0.035);
   h0->GetXaxis()->SetTitleSize(0.035);
   h0->GetXaxis()->SetTitleFont(42);
   h0->GetYaxis()->SetLabelFont(42);
   h0->GetYaxis()->SetLabelSize(0.035);
   h0->GetYaxis()->SetTitleSize(0.035);
   h0->GetYaxis()->SetTitleFont(42);
   h0->GetZaxis()->SetLabelFont(42);
   h0->GetZaxis()->SetLabelSize(0.035);
   h0->GetZaxis()->SetTitleSize(0.035);
   h0->GetZaxis()->SetTitleFont(42);
   h0->Draw("Same");
   
   TH1F *h1 = new TH1F("h1","ScintE {PID==1}",100,-400,4250);
   h1->SetBinContent(8,7);
   h1->SetBinContent(9,41594);
   h1->SetBinContent(10,4455);
   h1->SetBinContent(11,5062);
   h1->SetBinContent(12,4913);
   h1->SetBinContent(13,4488);
   h1->SetBinContent(14,4008);
   h1->SetBinContent(15,3277);
   h1->SetBinContent(16,2816);
   h1->SetBinContent(17,2208);
   h1->SetBinContent(18,1714);
   h1->SetBinContent(19,1472);
   h1->SetBinContent(20,1055);
   h1->SetBinContent(21,895);
   h1->SetBinContent(22,741);
   h1->SetBinContent(23,655);
   h1->SetBinContent(24,493);
   h1->SetBinContent(25,415);
   h1->SetBinContent(26,369);
   h1->SetBinContent(27,304);
   h1->SetBinContent(28,262);
   h1->SetBinContent(29,227);
   h1->SetBinContent(30,214);
   h1->SetBinContent(31,183);
   h1->SetBinContent(32,135);
   h1->SetBinContent(33,127);
   h1->SetBinContent(34,106);
   h1->SetBinContent(35,62);
   h1->SetBinContent(36,60);
   h1->SetBinContent(37,50);
   h1->SetBinContent(38,44);
   h1->SetBinContent(39,30);
   h1->SetBinContent(40,28);
   h1->SetBinContent(41,20);
   h1->SetBinContent(42,23);
   h1->SetBinContent(43,17);
   h1->SetBinContent(44,16);
   h1->SetBinContent(45,9);
   h1->SetBinContent(46,9);
   h1->SetBinContent(47,16);
   h1->SetBinContent(48,9);
   h1->SetBinContent(49,9);
   h1->SetBinContent(50,10);
   h1->SetBinContent(51,11);
   h1->SetBinContent(52,6);
   h1->SetBinContent(53,6);
   h1->SetBinContent(54,7);
   h1->SetBinContent(55,8);
   h1->SetBinContent(56,4);
   h1->SetBinContent(57,5);
   h1->SetBinContent(58,6);
   h1->SetBinContent(59,6);
   h1->SetBinContent(60,6);
   h1->SetBinContent(61,7);
   h1->SetBinContent(62,4);
   h1->SetBinContent(63,3);
   h1->SetBinContent(64,3);
   h1->SetBinContent(65,2);
   h1->SetBinContent(66,8);
   h1->SetBinContent(67,5);
   h1->SetBinContent(68,3);
   h1->SetBinContent(69,3);
   h1->SetBinContent(70,3);
   h1->SetBinContent(71,3);
   h1->SetBinContent(72,1);
   h1->SetBinContent(73,5);
   h1->SetBinContent(74,3);
   h1->SetBinContent(76,2);
   h1->SetBinContent(77,2);
   h1->SetBinContent(78,1);
   h1->SetBinContent(79,2);
   h1->SetBinContent(80,2);
   h1->SetBinContent(81,3);
   h1->SetBinContent(82,1);
   h1->SetBinContent(83,2);
   h1->SetBinContent(84,3);
   h1->SetBinContent(85,2);
   h1->SetBinContent(87,1);
   h1->SetBinContent(88,3);
   h1->SetBinContent(89,1);
   h1->SetBinContent(90,1);
   h1->SetBinContent(91,1);
   h1->SetBinContent(93,81);
   h1->SetEntries(82833);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   text = ptstats->AddText("h1");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 82833  ");
   text = ptstats->AddText("Mean  =  152.4");
   text = ptstats->AddText("RMS   =  275.8");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   h1->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(h1);
   h1->SetLineColor(2);
   h1->GetXaxis()->SetRange(1,100);
   h1->GetXaxis()->SetLabelFont(42);
   h1->GetXaxis()->SetLabelSize(0.035);
   h1->GetXaxis()->SetTitleSize(0.035);
   h1->GetXaxis()->SetTitleFont(42);
   h1->GetYaxis()->SetLabelFont(42);
   h1->GetYaxis()->SetLabelSize(0.035);
   h1->GetYaxis()->SetTitleSize(0.035);
   h1->GetYaxis()->SetTitleFont(42);
   h1->GetZaxis()->SetLabelFont(42);
   h1->GetZaxis()->SetLabelSize(0.035);
   h1->GetZaxis()->SetTitleSize(0.035);
   h1->GetZaxis()->SetTitleFont(42);
   h1->Draw("Same");
   
   TH1F *h2 = new TH1F("h2","ScintE {PID==2}",100,-400,4250);
   h2->SetBinContent(8,27);
   h2->SetBinContent(9,6528);
   h2->SetBinContent(10,479);
   h2->SetBinContent(11,262);
   h2->SetBinContent(12,197);
   h2->SetBinContent(13,150);
   h2->SetBinContent(14,147);
   h2->SetBinContent(15,129);
   h2->SetBinContent(16,119);
   h2->SetBinContent(17,172);
   h2->SetBinContent(18,214);
   h2->SetBinContent(19,202);
   h2->SetBinContent(20,189);
   h2->SetBinContent(21,229);
   h2->SetBinContent(22,221);
   h2->SetBinContent(23,195);
   h2->SetBinContent(24,179);
   h2->SetBinContent(25,157);
   h2->SetBinContent(26,161);
   h2->SetBinContent(27,145);
   h2->SetBinContent(28,146);
   h2->SetBinContent(29,139);
   h2->SetBinContent(30,126);
   h2->SetBinContent(31,115);
   h2->SetBinContent(32,103);
   h2->SetBinContent(33,98);
   h2->SetBinContent(34,101);
   h2->SetBinContent(35,91);
   h2->SetBinContent(36,89);
   h2->SetBinContent(37,111);
   h2->SetBinContent(38,91);
   h2->SetBinContent(39,89);
   h2->SetBinContent(40,81);
   h2->SetBinContent(41,72);
   h2->SetBinContent(42,71);
   h2->SetBinContent(43,81);
   h2->SetBinContent(44,65);
   h2->SetBinContent(45,63);
   h2->SetBinContent(46,62);
   h2->SetBinContent(47,51);
   h2->SetBinContent(48,51);
   h2->SetBinContent(49,48);
   h2->SetBinContent(50,36);
   h2->SetBinContent(51,57);
   h2->SetBinContent(52,39);
   h2->SetBinContent(53,31);
   h2->SetBinContent(54,33);
   h2->SetBinContent(55,46);
   h2->SetBinContent(56,30);
   h2->SetBinContent(57,38);
   h2->SetBinContent(58,33);
   h2->SetBinContent(59,31);
   h2->SetBinContent(60,31);
   h2->SetBinContent(61,33);
   h2->SetBinContent(62,25);
   h2->SetBinContent(63,23);
   h2->SetBinContent(64,23);
   h2->SetBinContent(65,22);
   h2->SetBinContent(66,18);
   h2->SetBinContent(67,20);
   h2->SetBinContent(68,17);
   h2->SetBinContent(69,17);
   h2->SetBinContent(70,12);
   h2->SetBinContent(71,18);
   h2->SetBinContent(72,15);
   h2->SetBinContent(73,19);
   h2->SetBinContent(74,16);
   h2->SetBinContent(75,10);
   h2->SetBinContent(76,11);
   h2->SetBinContent(77,14);
   h2->SetBinContent(78,13);
   h2->SetBinContent(79,9);
   h2->SetBinContent(80,7);
   h2->SetBinContent(81,11);
   h2->SetBinContent(82,11);
   h2->SetBinContent(83,13);
   h2->SetBinContent(84,12);
   h2->SetBinContent(85,15);
   h2->SetBinContent(86,9);
   h2->SetBinContent(87,9);
   h2->SetBinContent(88,8);
   h2->SetBinContent(89,15);
   h2->SetBinContent(90,8);
   h2->SetBinContent(91,3);
   h2->SetBinContent(93,413);
   h2->SetEntries(13320);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   text = ptstats->AddText("h2");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 13320  ");
   text = ptstats->AddText("Mean  =  579.4");
   text = ptstats->AddText("RMS   =  943.7");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   h2->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(h2);
   h2->SetLineColor(4);
   h2->GetXaxis()->SetRange(1,100);
   h2->GetXaxis()->SetLabelFont(42);
   h2->GetXaxis()->SetLabelSize(0.035);
   h2->GetXaxis()->SetTitleSize(0.035);
   h2->GetXaxis()->SetTitleFont(42);
   h2->GetYaxis()->SetLabelFont(42);
   h2->GetYaxis()->SetLabelSize(0.035);
   h2->GetYaxis()->SetTitleSize(0.035);
   h2->GetYaxis()->SetTitleFont(42);
   h2->GetZaxis()->SetLabelFont(42);
   h2->GetZaxis()->SetLabelSize(0.035);
   h2->GetZaxis()->SetTitleSize(0.035);
   h2->GetZaxis()->SetTitleFont(42);
   h2->Draw("Same");
   
   TH1F *h4 = new TH1F("h4","ScintE {PID==4}",100,-400,4250);
   h4->SetBinContent(8,99);
   h4->SetBinContent(9,35400);
   h4->SetBinContent(10,1536);
   h4->SetBinContent(11,127);
   h4->SetBinContent(12,56);
   h4->SetBinContent(13,50);
   h4->SetBinContent(15,2);
   h4->SetBinContent(24,1);
   h4->SetBinContent(28,1);
   h4->SetBinContent(33,1);
   h4->SetBinContent(41,23);
   h4->SetBinContent(42,136);
   h4->SetBinContent(43,143);
   h4->SetBinContent(44,115);
   h4->SetBinContent(45,124);
   h4->SetBinContent(46,102);
   h4->SetBinContent(47,121);
   h4->SetBinContent(48,76);
   h4->SetBinContent(49,79);
   h4->SetBinContent(50,81);
   h4->SetBinContent(51,73);
   h4->SetBinContent(52,62);
   h4->SetBinContent(53,64);
   h4->SetBinContent(54,70);
   h4->SetBinContent(55,55);
   h4->SetBinContent(56,72);
   h4->SetBinContent(57,63);
   h4->SetBinContent(58,57);
   h4->SetBinContent(59,71);
   h4->SetBinContent(60,55);
   h4->SetBinContent(61,77);
   h4->SetBinContent(62,94);
   h4->SetBinContent(63,90);
   h4->SetBinContent(64,110);
   h4->SetBinContent(65,124);
   h4->SetBinContent(66,146);
   h4->SetBinContent(67,187);
   h4->SetBinContent(68,242);
   h4->SetBinContent(69,367);
   h4->SetBinContent(70,458);
   h4->SetBinContent(71,552);
   h4->SetBinContent(72,704);
   h4->SetBinContent(73,856);
   h4->SetBinContent(74,964);
   h4->SetBinContent(75,1087);
   h4->SetBinContent(76,1145);
   h4->SetBinContent(77,1199);
   h4->SetBinContent(78,1093);
   h4->SetBinContent(79,1080);
   h4->SetBinContent(80,975);
   h4->SetBinContent(81,820);
   h4->SetBinContent(82,705);
   h4->SetBinContent(83,593);
   h4->SetBinContent(84,489);
   h4->SetBinContent(85,389);
   h4->SetBinContent(86,289);
   h4->SetBinContent(87,219);
   h4->SetBinContent(88,182);
   h4->SetBinContent(89,102);
   h4->SetBinContent(90,69);
   h4->SetBinContent(91,55);
   h4->SetBinContent(92,20);
   h4->SetBinContent(93,2143);
   h4->SetEntries(56540);
   
   ptstats = new TPaveStats(0.78,0.775,0.98,0.935,"brNDC");
   ptstats->SetName("stats");
   ptstats->SetBorderSize(1);
   ptstats->SetFillColor(0);
   ptstats->SetTextAlign(12);
   ptstats->SetTextFont(42);
   text = ptstats->AddText("h4");
   text->SetTextSize(0.0368);
   text = ptstats->AddText("Entries = 56540  ");
   text = ptstats->AddText("Mean  =   1064");
   text = ptstats->AddText("RMS   =   1512");
   ptstats->SetOptStat(1111);
   ptstats->SetOptFit(0);
   ptstats->Draw();
   h4->GetListOfFunctions()->Add(ptstats);
   ptstats->SetParent(h4);
   h4->SetLineColor(8);
   h4->GetXaxis()->SetRange(1,100);
   h4->GetXaxis()->SetLabelFont(42);
   h4->GetXaxis()->SetLabelSize(0.035);
   h4->GetXaxis()->SetTitleSize(0.035);
   h4->GetXaxis()->SetTitleFont(42);
   h4->GetYaxis()->SetLabelFont(42);
   h4->GetYaxis()->SetLabelSize(0.035);
   h4->GetYaxis()->SetTitleSize(0.035);
   h4->GetYaxis()->SetTitleFont(42);
   h4->GetZaxis()->SetLabelFont(42);
   h4->GetZaxis()->SetLabelSize(0.035);
   h4->GetZaxis()->SetTitleSize(0.035);
   h4->GetZaxis()->SetTitleFont(42);
   h4->Draw("Same");
   c1->Modified();
   c1->cd();
   c1->SetSelected(c1);
}
