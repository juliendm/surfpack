#! /usr/local/bin/python
from Tkinter import *
import math
import string
import os
import firstswig

root = Tk()

#first, a row for function entry and action button
fram = Frame(root)
Label(fram,text='dimensions').pack(side=LEFT)
ndims = Entry(fram, width = 3)
ndims.insert(END, "2")
ndims.pack(side=LEFT, fill=BOTH, expand=1)
ndim_button = Button(fram, text='Set')
ndim_button.pack(side=LEFT)
Label(fram,text='f(x)').pack(side=LEFT)
func = Entry(fram)
func.pack(side=LEFT, fill=BOTH, expand=1)
plot_button = Button(fram, text='Plot')
plot_button.pack(side=RIGHT)
fram.pack(side=TOP)

#then a row to enter bounds for f(x) in
#bounds = []
#for label in 'minX','maxX','minY','maxY':
#    Label(fram,text=label+':').pack(side=LEFT)
#    edit = Entry(fram, width = 6)
#    edit.pack(side=LEFT)
#    bounds.append(edit)


#Now add a row to hold the min and max values for f(x) in the plot
fram = Frame(root)
Label(fram,text="min f(x):").pack(side=LEFT)
minfxWidget = Entry(fram, width = 6)
minfxWidget.pack(side=LEFT)
minfxVal = -1.0
minfxWidget.insert(END, minfxVal)

Label(fram,text="max f(x):").pack(side=LEFT)
maxfxWidget = Entry(fram, width = 6)
maxfxWidget.pack(side=LEFT)
maxfxVal = 1.0
maxfxWidget.insert(END, maxfxVal)
fram.pack(side=TOP)
fram.pack_forget()
fram.pack(side=TOP)


# Number of dimensions in the data that is to be projected for plotting
currentNumDims = 2

# Array for constant value or min/max along dimensions that
# can be passed along the Python/C boundary. Entry pyPointDef[2n] is the
# min value for the nth dimension; pyPointDef[2n+1] is the max value for
# the nth dimension.
pyPointDefLength = currentNumDims * 2 
pyPointDef = firstswig.doubleArray(pyPointDefLength)
pyPointDef[0] = -1.0 	# default min x1
pyPointDef[1] = 1.0	# default max x1
pyPointDef[2] = -1.0	# default min x2
pyPointDef[3] = 1.0	# default max x2

# Add the widget to hold the fixed values/range for predictor variables
widgetForMinsAndMaxes = Frame(root)
widgetForMinsAndMaxes.pack(side = TOP)
activeDimVar = IntVar()
activeDimVar.set(0)

# Initialize structures for min/max widgets and values
dimframes = []
mlabels = []
minvals = []
maxvals = []
radiobuttons = []

def fillValues():
  global activeDimVar, maxvals, minvals, pyPointDef, pyPointDefLength
  for i in range(len(minvals)):
    minvals[i].delete(0,END)
    maxvals[i].delete(0,END)
    if i * 2 + 1 < pyPointDefLength:
      minvals[i].insert(END, pyPointDef[i*2])
      if activeDimVar.get() == i: 
        maxvals[i].insert(END, pyPointDef[i*2+1])
    else:
      print "Dim mismatch: minvals and pyPointDefLength"
      
    #  minvals[i].insert(END, -1.0)
    #  if activeDimVar.get() == i: 
    #    maxvals[i].insert(END, 1.0)

  
def radioSelection():
  global activeDimVar, mlabels, maxvals
  readValues()
  for i in range(currentNumDims):
    if activeDimVar.get() == i:
      mlabels[i].config(text="Range:")
      maxvals[i].config(borderwidth=2, highlightthickness=1,
        takefocus=True, state=NORMAL)
      if maxvals[i].cget('state') == 'normal': print "normal : ", i
      if i*2+1 < pyPointDefLength:
        maxvals[i].delete(0,END)
        maxvals[i].insert(END,pyPointDef[i*2+1])
    else:
      mlabels[i].config(text="Value:")
      maxvals[i].delete(0,END)
      maxvals[i].config(borderwidth=0,highlightthickness=0,
        takefocus=False, state=DISABLED)
      if maxvals[i].cget('state') == 'disabled': print "disabled : ", i

def readValues():
  print "Begin readValues"
  global pyPointDef, pyPointDefLength, minvals, maxvals
  #pyPointDefLength = len(minvals)*2
  #pyPointDef = firstswig.doubleArray(pyPointDefLength)
  for i in range(len(minvals)):
    try: pyPointDef[i*2] = float(minvals[i].get()) 
    except: 
      pyPointDef[i*2] = 0.0
      print "Bogus value in pyPointDef", i*2
    # Only attempt to read the max value if it is active
    if maxvals[i].cget('state') == 'normal': 
    #if activeDimVar.get() == i:
      try: pyPointDef[i*2+1] = float(maxvals[i].get()) 
      except: 
        print "Bogus value in pyPointDef", i*2+1, "|", maxvals[i].get(), "|", maxvals[i].cget('state')
        pyPointDef[i*2+1] = 0.0
  print "End readValues"

def removeFrames():
  global dimframes, mlabels, minvals, maxvals, radiobuttons
  #readValues()
  for ml in mlabels: ml.destroy()
  for mi in minvals: mi.destroy()
  for ma in maxvals: ma.destroy()
  for ra in radiobuttons: ra.destroy() 
  for l in dimframes: l.destroy()
  dimframes = []
  mlabels = []
  minvals = []
  maxvals = []
  radiobuttons = []
  

def setDimensionAction():
    print "Begin setDimensionAction"
    global currentNumDims, mlabels, pyPointDefLength, pyPointDef
    readValues()
    # store info about old values
    previousNumDims = currentNumDims
    previousPyPointDefLength = pyPointDefLength
    previousPyPointDef = pyPointDef
    for i in range(pyPointDefLength): print pyPointDef[i],
    print "Before"
    print "Previous number of dimensions: ", previousNumDims
    # acquire new number of dimensions from screen widget
    currentNumDims = int(ndims.get())
    print "New number of dimensions: ", currentNumDims
    if currentNumDims < 0: currentNumDims = 0
    print "New number of dimensions: ", currentNumDims
    # make a new array for the values
    pyPointDefLength = currentNumDims * 2
    pyPointDef = firstswig.doubleArray(pyPointDefLength)
    # copy values from before that are still pertinent
    for i in range(min(previousPyPointDefLength,pyPointDefLength)):
      pyPointDef[i] = previousPyPointDef[i]
    # If the number of dimensions increased, initialize new min/max values to (-1,1)
    for j in range(min(previousPyPointDefLength,pyPointDefLength),
		   max(previousPyPointDefLength,pyPointDefLength)):
      if j % 2 == 0:  # this is a min value
        pyPointDef[j] = -1.0
      else: # it's a max value
        pyPointDef[j] = 1.0 

    # If the number of dimensions decreased and the active dimension was removed,
    # update activeDimVar to first dimension
    for i in range(pyPointDefLength): print pyPointDef[i],
    print "After"
    if previousNumDims > currentNumDims: activeDimVar.set(0)
    removeFrames()
    addFrames()
    fillValues()
    radioSelection()
    print "End setDimensionAction"


def addFrames():
    for i in range(currentNumDims):
      dimframes.append( Frame(widgetForMinsAndMaxes) )
      mlabels.append( Label(dimframes[i], text = 'Value: '))
      mlabels[i].pack(side=LEFT)
      minvals.append( Entry(dimframes[i], width=5) )
      minvals[i].pack(side=LEFT)
      maxvals.append( Entry(dimframes[i], width=5))
      maxvals[i].pack(side=LEFT)
      if activeDimVar.get() != i:
        maxvals[i].config(state = DISABLED)
      radiobuttons.append( Radiobutton(dimframes[i],variable=activeDimVar,value=i,
	 command=radioSelection))
      radiobuttons[i].pack()
      dimframes[i].pack()
    #activeDimVar.set(0)
    #fillValues()
    #radioSelection()

#execute initialization
#setDimensionAction()
addFrames()
fillValues()
radioSelection()
    
    

#and finally the canvas
c = Canvas(root)
c.pack(side=TOP, fill=BOTH, expand=1)

cedit = Entry(root)
cedit.pack(side=TOP, fill=BOTH, expand=1)

myobj = firstswig.FirstClass(4)
k = 1000
def enterPressedCommand(*ignore):
    global myobj
    global k
    str = cedit.get()
    print str
    t.insert(END, str + "\n")
    cedit.delete(0,END)
    sed_command = "cat template | sed 's/\/\/REPLACE_ME/" + str \
        + "/' > firstswig.cxx"
    print sed_command
    os.system(sed_command)
    os.system("./mymake")
    #import firstswig
    #reload(firstswig._firstswig)
    reload(firstswig)
    
    myobj = firstswig.FirstClass(k)
    k = k*2
    myobj.printVal("Will it work?")
    
    myarray =  firstswig.doubleArray(3)
    myarray[0] = 1.0
    myarray[1] = 2.0
    myarray[2] = 4.0
    
    #x = myobj.shiftArray(myarray, 3)
    #print x
    #x = myobj.shiftArray(myarray, 3)
    #print x
    y= myobj.myevaluate(myarray, 3)
    print y

    

cedit.bind('<Return>',enterPressedCommand)

t = Text(root)
t.pack(side=TOP, fill=BOTH, expand=1)
t.start_command = '1.0'
def enterPressed(*ignore): 
    s = t.get(t.start_command, END)
    print "[", string.strip(s), "]"
    print t.start_command
    t.insert(END, "\n_____ ")
    t.start_command = CURRENT

t.bind('<Return>', enterPressed)

#def minimax(values=[0.0, 1.0, 0.0, 1.0]):
#    "Adjust and display X and Y bounds"
#    for i in range(4):
#        edit = bounds[i]
#        try: values[i] = float(edit.get())
#        except: pass
#        edit.delete(0, END)
#        edit.insert(END, '%.2f'%values[i])
#    return values

def plot():
    "Plot given functions with given bounds"
    #minx, maxx, miny, maxy = minimax()
    minx, maxx, miny, maxy =[0.0, 10.0, -2.0, 2.0]
    
    #get and compile the function
    f = func.get()
    f = compile(f, f, 'eval')
    
    #get Canvas X and Y dimensions
    CX = c.winfo_width()
    CY = c.winfo_height()

    #compute coordinates for line
    coords = []
    for i in range(0,CX,5):
        coords.append(i)
        x = minx + ((maxx-minx)*i)/CX
        z = 5
        y = eval(f, vars(math), {'x':x,'z':z})
        j = CY - CY*(y-miny)/(maxy - miny)
        coords.append(j)

    # draw line
    c.delete(ALL)
    c.create_line(*coords)
 
plot_button.config(command=plot)
ndim_button.config(command=setDimensionAction)
#minvals[0].config(state=DISABLED)


# give an initial example in lieu of docs
# Initial value for the function text field
# This is the function that will be plotted.
f = 'sin(x) + cos(x)'
func.insert(END, f)
#minimax([0.0, 10.0, -2.0, 2.0])

def gotit(var) :
  print var
  return var + 3

# Initialize command field with an example command
cedit.delete(0,END)
cedit.insert(END, "result = vals[0] + vals[1];")
cedit.focus_set()
if __name__ == '__main__' :
  s = input("X: ")
  gotit(s)
  try: root.mainloop()
  except: print "Error"
