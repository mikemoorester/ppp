import unittest 

import residuals as res

class parseDPH(unittest.TestCase):

    # will it parse these files
    def testOne(self):
        self.failUnless(res.parseDPH('./t/DPH.TOW2.001.PRNAL.gz'))
        self.failUnless(res.parseDPH('./t/DPH.TOW2.001.PRNAL'))
        self.failUnless(res.parseDPH('./t/DPH.TOW2.004.PRNAL.gz'))

    def testFileParse(self):
        #dphs1 = res.parseDPH('./t/DPH.TOW2.001.PRNAL.gz')
        dphs1 = res.parseDPH('./t/DPH.TOW2.001.PRNAL')
        # test the correct keys exists
        self.failUnless('satsViewed' in dphs1)
        self.failUnless('epochs' in dphs1)
        
        # check that each satellite viewed has a key in the dictionary
        self.failUnless(31 == len(dphs1['satsViewed']))
        for sat in dphs1['satsViewed']:
            prn = 'prn_'+str(sat)
            self.failUnless(prn in dphs1)

        # check that each epoch observed has a key
        self.failUnless(2880 == len(dphs1['epochs']))

        for epoch in dphs1['epochs']:
            self.failUnless(str(epoch) in dphs1)

        # Test getting the data from one epoch
        self.failUnless(13 == len(dphs1['1']))

        # check that each observation has the following keys defined
        satObsKeys = ['l1cyc', 'pccyc', 'p2cyc', 'ncyc', 'l2cyc', 'dataf', 'prn',
                'p1cyc', 'lsv', 'epoch', 'lgcyc','pf','wlcyc', 'lccyc', 'el', 'az']

        print(dphs1['781'])
        print(dphs1['781'][4])

        self.failUnless(1 in dphs1['epochs'])

        for epoch in dphs1['epochs']:
            print('epoch:',epoch)
            print('dphs1[epochStr]:',dphs1[str(epoch)])
            #ep = dphs1[str(epoch)][sat] #int(epoch) - 1
            #ep = int(epoch)
            #epoch = ep
            #for sat in dphs1['1']:
            for sat in dphs1[str(epoch)]:
                #print('Sat:',sat)
                ep = dphs1[str(epoch)][str(sat)]
                print('ep',ep)
                satkey = 'prn_'+str(sat)
                for key in satObsKeys:
                    self.failUnless(key in dphs1[satkey][ep])

def main():
    unittest.main()

if __name__ == '__main__':
    main()
