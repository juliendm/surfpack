/*  _______________________________________________________________________
 
    Surfpack: A Software Library of Multidimensional Surface Fitting Methods
    Copyright (c) 2006, Sandia National Laboratories.
    This software is distributed under the GNU Lesser General Public License.
    For more information, see the README file in the top Surfpack directory.
    _______________________________________________________________________ */

//
// "$Id: spgui.cpp 454 2009-12-18 23:57:00Z briadam $"
//
// OpenGL overlay test program for the Fast Light Tool Kit (FLTK).
//
// Copyright 1998-2005 by Bill Spitzak and others.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public
// License as published by the Free Software Foundation; either
// version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA.
//
// Please report all bugs and problems on the following page:
//
//     http://www.fltk.org/str.php
//

#include "config.h"
// Surfpack headers
#include "surfpack.h"
#include "SurfpackInterpreter.h"
#include "SurfpackParser.h"
#include "SurfpackParserArgs.h"
#include "Surface.h"
#include "SurfData.h"
#include "SurfPoint.h"
#include "SurfaceFactory.h"

// std C/C++ headers
#include <stdio.h>
#include <vector>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <cassert>

// FLTK headers
#include <FL/fl_ask.H>
#include <FL/Fl.H>
#include <FL/Fl_Window.H>
#include <FL/Fl_Input.H>
#include <FL/Fl_Int_Input.H>
#include <FL/Fl_Float_Input.H>
#include <FL/Fl_Output.H>
#include <FL/Fl_Button.H>
#include <FL/Fl_Hor_Slider.H>
#include <FL/Fl_Value_Slider.H>
#include <FL/Fl_Toggle_Button.H>
#include <FL/Fl_Menu_Bar.H>
#include <FL/Fl_Box.H>
#include <FL/math.h>


// Global Data
typedef std::vector<double> Pt;
typedef std::vector<double> ResponseList;
Pt pta;
Pt ptb;
std::vector<int> actives(4);
std::vector<Surface*> surfs(4);
SurfData* g_data;
ResponseList rl(4);
float colors[4][3] =  {{1,0,0},{0,1,0},{0,0,1},{1,1,1}};

void currPt(Pt& pt, const Pt& x, const Pt& y, double alpha)
{
  assert(alpha >= 0.0 && alpha <= 1.001);
  assert(x.size() == y.size());
  pt.resize(x.size());
  for (unsigned i = 0; i < pt.size(); i++) {
    pt[i] = alpha*y[i]+(1.0-alpha)*x[i];
  }
}

class Evaluator {
public:
  Evaluator() {}
  virtual double eval(double* pt, unsigned size) = 0;
};
class SurfaceEvaluator : public Evaluator
{
public:
  SurfaceEvaluator() : s(0) {}
private:
  Surface* s;
public:
  void set(Surface* s_) { s = s_;}
  double eval(double* pt, unsigned size) {
    assert(s);
    std::vector<double> spt;
    for (unsigned i = 0; i < size ; i++) {
      spt.push_back(pt[i]);
    }
    return s->getValue(spt);
  } 
};

class Sine : public Evaluator
{
public:
  Sine() {}
  double eval(double* pt, unsigned size)
  {
    double res  = 2;
    for (unsigned i = 0; i < size; i++) {
      res += sin(pt[i]);
    }
    return res;
  }
};
    

class Cosine : public Evaluator
{
public:
  Cosine() {}
  double eval(double* pt, unsigned size)
  {
    double res  = 0;
    for (unsigned i = 0; i < size; i++) {
      res += cos(pt[i]);
    }
    return res;
  }
};

class Tanh : public Evaluator
{
public:
  Tanh() {}
  double eval(double* pt, unsigned size)
  {
    double res  = 0;
    for (unsigned i = 0; i < size; i++) {
      res += tanh(pt[i]);
    }
    return res;
  }
};

class Zero : public Evaluator
{
public:
  Zero() {}
  double eval(double* pt, unsigned size)
  {
    return 0.0;
  }
};

std::vector<Evaluator*> funcs(4);

void init()
{
  //funcs[0] = new Sine; 
  //funcs[1] = new Cosine; 
  //funcs[2] = new Tanh; 
  //funcs[3] = new Zero; 
  for (unsigned i = 0; i < 4; i++) {
    surfs[i] = 0;
    funcs[i] = new SurfaceEvaluator();
  }
}

void cleanup()
{
  for (unsigned i = 0; i < funcs.size(); i++) {
    delete funcs[i];
    delete surfs[i];
  }
}

void parsePt(Pt& pt, const char* charpt)
{
  pt.clear();
  std::string s(charpt);
  std::istringstream is(s);
  double element;
  while(true) {
    is >> element;
    if (!is.fail()) {
      pt.push_back(element);
    }
    if (is.eof()) break;
  }
  copy(pt.begin(),pt.end(),std::ostream_iterator<double>(std::cout," "));
  std::cout << std::endl;
}

#if !HAVE_GL
#include <FL/Fl_Box.H>
class shape_window : public Fl_Box {
public:	
  int sides;
  shape_window(int x,int y,int w,int h,const char *l=0)
    :Fl_Box(FL_DOWN_BOX,x,y,w,h,l){
      label("This demo does\nnot work without GL");
  }
};
#else
#include <FL/gl.h>
#include <FL/Fl_Gl_Window.H>

// global variables
float min;
float max;
int intervals;
std::vector<ResponseList> evals;
Fl_Float_Input* g_range_min;
Fl_Float_Input* g_range_max;
Fl_Toggle_Button* g_auto_range;

void do_evaluations()
{
  if (intervals <= 0) return; 
  double minf = 1e10;
  double maxf = -1e10;
  double inc = 1.0 / intervals;
  evals.clear();
  double alpha = 0.0;
  double curf;
  Pt pt;
  assert(funcs.size() == actives.size());
  do {
    //printf("alpha: %f\n",alpha);
    currPt(pt,pta,ptb,alpha);
    copy(pt.begin(),pt.end(),std::ostream_iterator<double>(std::cout,","));
    for (unsigned i = 0; i < funcs.size(); i++) {
      if (actives[i]) {
        curf = rl[i] = funcs[i]->eval(&pt[0],pt.size());
        std::cout << rl[i] << "|" ;
        if (curf > maxf) maxf = curf;
        if (curf < minf) minf = curf;
      }
    }
    std::cout << std::endl;
    evals.push_back(rl);
    alpha += inc;
  } while (evals.size() <= (unsigned)intervals);
  if (g_auto_range->value()) {
    std::ostringstream os; os << minf;
    g_range_min->value(os.str().c_str());
    std::ostringstream os2; os2 << maxf;
    g_range_max->value(os2.str().c_str());
  }
}


class shape_window : public Fl_Gl_Window {
  void draw();
  void draw_overlay();
public:
  int sides;
  int overlay_sides;
  shape_window(int x,int y,int w,int h,const char *l=0);
};

shape_window::shape_window(int x,int y,int w,int h,const char *l) :
Fl_Gl_Window(x,y,w,h,l) {
  sides = overlay_sides = 3;
}

void shape_window::draw() {
// the valid() property may be used to avoid reinitializing your
// GL transformation for each redraw:
  if (!valid()) {
    valid(1);
    glLoadIdentity();
    glViewport(0,0,w(),h());
  }
  double minf = atof(g_range_min->value());
  double maxf = atof(g_range_max->value());
// draw an amazing but slow graphic:
  glClear(GL_COLOR_BUFFER_BIT);
  //  for (int j=1; j<=1000; j++) {
    glLineWidth(3.0);
    for (unsigned j = 0; j < actives.size(); j++) {
      if (actives[j] and surfs[j]) {
        glColor3f(colors[j][0],colors[j][1],colors[j][2]);
        for (unsigned i=0; i < evals.size(); i++) {
          float xval = (i/(float)evals.size())*2.0f - 1.0f;
          float yval = ((evals[i][j]-minf)/(maxf-minf))*2.0f - 1.0f;
          printf("%f %f\n",xval,yval);
        }
        glBegin(GL_LINE_STRIP);
        for (unsigned i=0; i< evals.size(); i++) {
          float xval = (i/(float)(evals.size()-1))*2.0f - 1.0f;
          float yval = ((evals[i][j]-minf)/(maxf-minf))*2.0f - 1.0f;
          //glColor3f(1.0,1.0,1.0);
          glVertex3f(xval,yval,0);
        }
        glEnd();
      }
    }
  // }
}

void shape_window::draw_overlay() {
// the valid() property may be used to avoid reinitializing your
// GL transformation for each redraw:
  return;
  if (!valid()) {
    valid(1);
    glLoadIdentity();
    glViewport(0,0,w(),h());
  }
// draw an amazing graphic:
  gl_color(FL_RED);
  glBegin(GL_LINE_LOOP);
  for (int i=0; i<overlay_sides; i++) {
    double ang = i*2*M_PI/overlay_sides;
    glVertex3f(cos(ang),sin(ang),0);
  }
  glEnd();
}
#endif

// when you change the data, as in this callback, you must call redraw():

// CALLBACKS 

void load_cb(Fl_Widget* o, void*)
{
  const char* filename = fl_input("Surfpack data filename","quad.spd");
  try {
    g_data = new SurfData(filename);
  } catch(...) {
    fl_message("Error opening data file.");
    return;
  } 
  try {
    for (unsigned i = 0; i < 4; i++) {
      delete surfs[i];
      switch (i) {
        case 0: 
          surfs[i] = SurfaceFactory::createSurface("ann",g_data);
          break;
        case 1: 
          surfs[i] = SurfaceFactory::createSurface("kriging",g_data);
          break;
        case 2: 
        default:
          surfs[i] = SurfaceFactory::createSurface("mars",g_data);
          break;
        case 3: 
          surfs[i] = SurfaceFactory::createSurface("polynomial",g_data);
          break;
      }
      surfs[i]->createModel();
      dynamic_cast<SurfaceEvaluator*>(funcs[i])->set(surfs[i]);
    } // for
  } catch (...) {
    fl_message("Error creating surfaces");
    // Possibly leak some memory
    for (unsigned i = 0; i < 4; i++) {
      surfs[i] = 0;
    }
  }
}
void auto_range_cb(Fl_Widget* o)
{
  if (((Fl_Toggle_Button*)o)->value()) {
    g_range_min->deactivate();
    g_range_max->deactivate();
  } else {
    g_range_min->activate();
    g_range_max->activate();
  }
}

void vslider_cb(Fl_Widget* o)
{
  intervals = (int)((Fl_Value_Slider*)o)->value();
  printf("intervals: %d\n",intervals);
}
void graph_cb(Fl_Widget* o, void* arg)
{
  do_evaluations();
  shape_window *sw = (shape_window *)arg;
  sw->redraw();
  printf("graph callback\n");
}

void menu_cb(Fl_Widget* o)
{
  Fl_Menu_* mw = (Fl_Menu_*)o;
  const Fl_Menu_Item* m = mw->mvalue();
  int user_data = (int)(m->user_data());
  actives[user_data] = m->value();
  printf("%d %s %d\n",user_data, m->label(),m->value()); 
}

void basic_cb(Fl_Widget* o)
{
  printf("Callback for %s '%s'\n",o->label(),((Fl_Input*)o)->value());
}

void sides_cb(Fl_Widget *o, void *p) {
  shape_window *sw = (shape_window *)p;
  sw->sides = int(((Fl_Slider *)o)->value());
  sw->redraw();
}

const char* global_string;
void button_cb(Fl_Widget* o, void* arg)
{
  printf("button_cb callback: %s\n",global_string);
  ((Fl_Output*)arg)->value(global_string);
}

void bpoint_cb(Fl_Widget* o)
{
  global_string = ((Fl_Input*)o)->value();
  //printf("global: %s\n",global_string);
  //sscanf(global_string,"%f %f %d",&min,&max,&intervals);
  //printf("min %f max %f intervals %d\n",min,max,intervals);
  //do_evaluations();
  parsePt(pta,((Fl_Input*)o)->value());
}

void epoint_cb(Fl_Widget* o)
{
  parsePt(ptb,((Fl_Input*)o)->value());
}

void copy_cb(Fl_Widget* o, void* arg)
{
  ((Fl_Input*)arg)->value(global_string);
}

#if HAVE_GL
void overlay_sides_cb(Fl_Widget *o, void *p) {
  shape_window *sw = (shape_window *)p;
  sw->overlay_sides = int(((Fl_Slider *)o)->value());
  sw->redraw_overlay();
}
#endif
#include <stdio.h>


Fl_Menu_Item menutable[] = {
  {"&Surface",0,0,0,FL_SUBMENU},
    {"&ann",FL_CTRL+'a',0,(void*)0,FL_MENU_TOGGLE},
    {"&kriging",FL_CTRL+'k',0,(void*)1,FL_MENU_TOGGLE},
    {"&mars",FL_CTRL+'m',0,(void*)2,FL_MENU_TOGGLE},
    {"&quadratic",FL_CTRL+'q',0,(void*)3,FL_MENU_TOGGLE},
    {0},
  {"&Data",0,0,0,FL_SUBMENU},
    {"&Load",FL_CTRL+'l',load_cb,0,0},
    {0},
  {0}
};

int main(int argc, char **argv) {
  intervals = 10;
  init();
  const int WIDTH=400;
  Fl_Window window(WIDTH, 500);

  shape_window sw(10, 165, window.w()-20, window.h()-180);
//sw.mode(FL_RGB);
  window.resizable(&sw);

  Fl_Value_Slider vslider(60, 40, window.w()-70, 30, "Points:");
  vslider.align(FL_ALIGN_LEFT);
  vslider.step(1);
  vslider.precision(0);
  vslider.type(FL_HORIZONTAL);
  vslider.bounds(1,500);
  vslider.value(10);
  vslider.callback(vslider_cb);

  // Add begin field
  Fl_Input begin_point(60, 75, window.w()-70,30, "Begin:");
  begin_point.tooltip("Enter the values, separated by spaces.");
  begin_point.callback(bpoint_cb);
  begin_point.value("0");
  bpoint_cb(&begin_point);

  // Output field
  Fl_Input end_point(60, 105, window.w()-70,30, "End: ");
  end_point.callback(epoint_cb);
  end_point.value("5");
  epoint_cb(&end_point);

  // Numpoints field
  Fl_Float_Input range_min(60, 135, 85, 30, "Range:");
  g_range_min = &range_min;
  Fl_Float_Input range_max(150, 135, 85, 30, 0);
  g_range_max = &range_max;

  // Copy Button 
  Fl_Toggle_Button auto_range(245, 135, 70, 30, "Auto");
  g_auto_range = &auto_range;
  auto_range.callback(auto_range_cb);
  auto_range.set();
  
  Fl_Button graph_button(320, 135, 70, 30, "Graph");
  graph_button.callback(graph_cb,&sw);
  Fl_Button copy_button(320, 135, 0, 0, 0);
  //copy_button.hide();
  copy_button.shortcut(FL_CTRL+'c');
  copy_button.callback(copy_cb,&end_point);

  Fl_Menu_Bar mb(0,0,WIDTH,30,0);
  mb.menu(menutable);
  mb.callback(menu_cb);

  window.end();
  window.show(argc,argv);
#if HAVE_GL
  printf("Can do overlay = %d\n", sw.can_do_overlay());
  sw.show();
  sw.redraw_overlay();
#else
  printf("Cannot do overlay\n");
  sw.show();
#endif

  return Fl::run();
  cleanup();
}

//
// End of "$Id: spgui.cpp 454 2009-12-18 23:57:00Z briadam $".
//
