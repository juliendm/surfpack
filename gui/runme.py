import firstswig

myobj = firstswig.FirstClass(4)
#myobj.printVal("Will it work?")
#
myarray =  firstswig.doubleArray(3)
myarray[0] = 1.0
myarray[1] = 2.0
myarray[2] = 4.0
#
x = myobj.shiftArray(myarray, 3)
print x
#x = myobj.shiftArray(myarray, 3)
#print x

y = myobj.myevaluate(myarray, 3)
print y

