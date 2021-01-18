#include "TStyle.h"

//// tdrGrid: Turns the grid lines on (true) or off (false)
//void tdrGrid(bool gridOn)
//{
//  tdrStyle->SetPadGridX(gridOn);
//  tdrStyle->SetPadGridY(gridOn);
//}
//
//// fixOverlay: Redraws the axis
//void fixOverlay()
//{
//  gPad->RedrawAxis();
//}

void setTDRStyle(bool show_statsbox=true)
{
    TStyle* tdrStyle = new TStyle("tdrStyle","Style for P-TDR");

    Int_t font_num = 42; // 132
    
    // For the Canvas:
    tdrStyle->SetCanvasBorderMode(1);
    tdrStyle->SetCanvasColor(kWhite);
    tdrStyle->SetCanvasDefH(620); //Height of canvas
    tdrStyle->SetCanvasDefW(600); //Width of canvas
    tdrStyle->SetCanvasDefX(0);   //POsition on screen
    tdrStyle->SetCanvasDefY(0);

    // For the Pad:
    tdrStyle->SetPadBorderMode(0);
//  tdrStyle->SetPadBorderSize(Width_t size = 1);
    tdrStyle->SetPadColor(kWhite);
    tdrStyle->SetPadGridX(false);
    tdrStyle->SetPadGridY(false);
    tdrStyle->SetGridColor(0);
    tdrStyle->SetGridStyle(3);
    tdrStyle->SetGridWidth(1);

    // For the frame:
    tdrStyle->SetFrameBorderMode(0);
    tdrStyle->SetFrameBorderSize(1);
    tdrStyle->SetFrameFillColor(0);
    tdrStyle->SetFrameFillStyle(0);
    tdrStyle->SetFrameLineColor(1);
    tdrStyle->SetFrameLineStyle(1);
    tdrStyle->SetFrameLineWidth(1);

    // For the histo:
    // tdrStyle->SetHistFillColor(1);  // 0
    // tdrStyle->SetHistFillStyle(0);
    tdrStyle->SetHistLineColor(4); // 1
    tdrStyle->SetHistLineStyle(0);
    tdrStyle->SetHistLineWidth(1);
//  tdrStyle->SetLegoInnerR(Float_t rad = 0.5);
    tdrStyle->SetNumberContours(Int_t number = 200);

    tdrStyle->SetEndErrorSize(2);
//  tdrStyle->SetErrorMarker(20);
    // tdrStyle->SetErrorX(0.);

    tdrStyle->SetMarkerStyle(20);

    //For the fit/function:
    tdrStyle->SetOptFit(1);
    tdrStyle->SetFitFormat("5.4g");
    tdrStyle->SetFuncColor(2);
    tdrStyle->SetFuncStyle(1);
    tdrStyle->SetFuncWidth(1);

    //For the date:
    tdrStyle->SetOptDate(0);
// tdrStyle->SetDateX(Float_t x = 0.01);
// tdrStyle->SetDateY(Float_t y = 0.01);

    // For the statistics box:
    tdrStyle->SetOptFile(0);
    if show_statsbox {
        tdrStyle->SetOptStat("iouRMe"); 
    } else {
        tdrStyle->SetOptStat(0);
        // For the fit function:
        tdrStyle->SetOptFit(0);
    }
    tdrStyle->SetStatColor(kWhite);
    tdrStyle->SetStatFont(font_num);
    tdrStyle->SetStatFontSize(0.025);  // 0.03
    tdrStyle->SetStatTextColor(1);
    tdrStyle->SetStatFormat(".4f");  // 6.4g width?.<num_sig_figs>
    tdrStyle->SetStatBorderSize(1);
    // tdrStyle->SetStatH(0.1);  // 1, 0.05 is pretty good. 0.02 is too small
    // tdrStyle->SetStatW(0.15);  // # 0.15 is pretty good.
    // tdrStyle->SetStatStyle(Style_t style = 1001);
    tdrStyle->SetStatX(Float_t x = 0.90);
    tdrStyle->SetStatY(Float_t y = 0.90);

    // Margins:
    tdrStyle->SetPadTopMargin(0.07);
    tdrStyle->SetPadBottomMargin(0.17);
    tdrStyle->SetPadLeftMargin(0.15);
    tdrStyle->SetPadRightMargin(0.05);

    // For the Global title:
    tdrStyle->SetOptTitle(1);  // 0 is "off": no title.
    // tdrStyle->SetTitleColor(1);
    // tdrStyle->SetTitleTextColor(1);
    tdrStyle->SetTitleFillColor(10);
    tdrStyle->SetTitleFontSize(0.031);  // 0.05
    tdrStyle->SetTitleFont(font_num, "somethingelse");  // 42, Anything other "XYZ" will change the title.
// tdrStyle->SetTitleH(0); // Set the height of the title box
    tdrStyle->SetTitleW(1.0); // Set the width of the title box, 0.8
    tdrStyle->SetTitleX(0.0); // Set the position of the title box, 0.1
    tdrStyle->SetTitleY(0.965); // Set the position of the title box, 0.985
// tdrStyle->SetTitleStyle(Style_t style = 1001);
    tdrStyle->SetTitleBorderSize(0);
    tdrStyle->SetTitleAlign(13);  // 11 is against left side, 22 is centered vert and horiz, 21

    // For the axis titles:
    tdrStyle->SetTitleColor(1, "XYZ");
    tdrStyle->SetTitleFont(font_num, "XYZ");  // 42
    tdrStyle->SetTitleSize(0.030, "XYZ");  // 0.06
// tdrStyle->SetTitleXSize(Float_t size = 0.02); // Another way to set the size?
// tdrStyle->SetTitleYSize(Float_t size = 0.02);
    tdrStyle->SetTitleXOffset(1.1);  // 1.0
    tdrStyle->SetTitleYOffset(1.3);  // 1.1
// tdrStyle->SetTitleOffset(1.1, "Y"); // Another way to set the Offset

    // For the axis labels:
    tdrStyle->SetLabelColor(1, "XYZ");
    tdrStyle->SetLabelFont(font_num, "XYZ");  // This is font style, not size!
    tdrStyle->SetLabelOffset(0.007, "XYZ");
    tdrStyle->SetLabelSize(0.02, "XYZ");  // 0.05

    // For the axis:
    tdrStyle->SetAxisColor(1, "XYZ");
    tdrStyle->SetStripDecimals(kTRUE);
    tdrStyle->SetTickLength(0.03, "XYZ");
    tdrStyle->SetNdivisions(510, "XYZ");
    tdrStyle->SetPadTickX(1);  // To get tick marks on the opposite side of the frame
    tdrStyle->SetPadTickY(1);

    // For the legend:
    tdrStyle->SetLegendFont(font_num);
    tdrStyle->SetLegendTextSize(0.025);

    // Change for log plots:
    tdrStyle->SetOptLogx(0);
    tdrStyle->SetOptLogy(0);
    tdrStyle->SetOptLogz(0);

    // Postscript options:
    tdrStyle->SetPaperSize(20.,20.);
// tdrStyle->SetLineScalePS(Float_t scale = 3);
// tdrStyle->SetLineStyleString(Int_t i, const char* text);
// tdrStyle->SetHeaderPS(const char* header);
// tdrStyle->SetTitlePS(const char* pstitle);

// tdrStyle->SetBarOffset(Float_t baroff = 0.5);
// tdrStyle->SetBarWidth(Float_t barwidth = 0.5);
// tdrStyle->SetPaintTextFormat(const char* format = "g");
// tdrStyle->SetPalette(Int_t ncolors = 0, Int_t* colors = 0);
// tdrStyle->SetTimeOffset(Double_t toffset);
// tdrStyle->SetHistMinimumZero(kTRUE);
    tdrStyle->SetHatchesLineWidth(5);
    tdrStyle->SetHatchesSpacing(0.05);

    tdrStyle->cd();
}
