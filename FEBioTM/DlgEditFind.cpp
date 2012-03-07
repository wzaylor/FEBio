#include "stdafx.h"
#include "DlgEditFind.h"
#include <FL/Fl_Input.H>
#include <FL/Fl_Check_Button.H>
#include <FL/Fl_Return_Button.H>
#include <Flex.h>

CDlgEditFind::CDlgEditFind() : Flx_Dialog(400, 150, "Find")
{
	int W = w();
	int H = h();

	m_sztxt[0] = 0;
	m_bcase = false;

	Fl_Input* pt;
	Fl_Check_Button* pc;
	Fl_Button* pb;
	Fl_Return_Button* pr;
	begin();
	{
		pt = new Fl_Input(5, 30, W-10, 20, "Find:"); pt->align(FL_ALIGN_LEFT | FL_ALIGN_TOP); DDX_TEXT(pt, m_sztxt);
		pc = new Fl_Check_Button(5, 60, 100, 20, "case sensitive"); DDX_BOOL(pc, &m_bcase);

		pr = new Fl_Return_Button(W/2-75, H-40, 70, 30, "Find"  ); AddCallback(pr, (FLX_CALLBACK) &Flx_Dialog::OnOk);
		pb = new Fl_Button       (W/2   , H-40, 70, 30, "Cancel"); AddCallback(pb, (FLX_CALLBACK) &Flx_Dialog::OnCancel);
	}
	end();
}

CDlgEditGoToLine::CDlgEditGoToLine() : Flx_Dialog(200, 100, "Go to line")
{
	int W = w();
	int H = h();

	m_nline = 1;

	Flx_Int_Input* pi;
	Fl_Button* pb;
	Fl_Return_Button* pr;
	begin();
	{
		pi = new Flx_Int_Input(50, 20, W-60, 20, "Line:"); pi->align(FL_ALIGN_LEFT); DDX_INT(pi, &m_nline);

		pr = new Fl_Return_Button(W/2-75, H-40, 70, 30, "OK"    ); AddCallback(pr, (FLX_CALLBACK) &Flx_Dialog::OnOk);
		pb = new Fl_Button       (W/2   , H-40, 70, 30, "Cancel"); AddCallback(pb, (FLX_CALLBACK) &Flx_Dialog::OnCancel);
	}
	end();
}
