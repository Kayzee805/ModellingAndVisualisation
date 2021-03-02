import test
import time

def main():
    size = 50
    sweeps = 10000
    initialisation="random"
    t1=time.time()
   # test.animate(size,sweeps,initialisation)
    test.generateHistogram(size,sweeps,initialisation)
  #  test.generateCom(size,sweeps,initialisation)

    t2=time.time()
    print(f"Time taken for all = {t2-t1}s")

    #test.plotAll()
main()