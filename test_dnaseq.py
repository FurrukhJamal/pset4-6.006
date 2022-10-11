from dnaseq import *

### Testing ###

class TestRollingHash(unittest.TestCase):
    def test_rolling(self):
        rh1 = RollingHash('CTAGC')
        rh2 = RollingHash('TAGCG')
        rh3 = RollingHash('AGCGT')
        rh1.slide('C','G')
        self.assertTrue(rh1.current_hash() == rh2.current_hash())
        rh1.slide('T','T')
        self.assertTrue(rh1.current_hash() == rh3.current_hash())

class TestMultidict(unittest.TestCase):
    def test_multi(self):
        foo = Multidict()
        foo.put(1, 'a')
        foo.put(2, 'b')
        foo.put(1, 'c')
        self.assertTrue(foo.get(1) == ['a','c'])
        self.assertTrue(foo.get(2) == ['b'])
        self.assertTrue(foo.get(3) == [])

# This test case may break once you add the argument m (skipping).
class TestExactSubmatches(unittest.TestCase):
   def test_one(self):
       foo = 'yabcabcabcz'
       bar = 'xxabcxxxx'
       matches = list(getExactSubmatches(iter(foo), iter(bar), 3, 1))
       correct = [(1,2), (4,2), (7,2)]
       self.assertTrue(len(matches) == len(correct))
       for x in correct:
           self.assertTrue(x in matches)


class TestintervalSubsequenceHashes(unittest.TestCase):
    def test_interval(self):
        foo = "abcxyzabcdefabcklmabcku"
        matches = list(intervalSubsequenceHashes(iter(foo), 3, 3))
        print(f"matches : {matches}")
        rhash = RollingHash("abc")
        hash = rhash.current_hash()
        correct = [(hash , ("abc", 0)), (hash , ("abc", 6)), (hash , ("abc", 12)), (hash , ("abc", 18))]
        self.assertTrue(len(matches) == len(correct))
        self.assertTrue(matches[0][0] == correct[0][0])
        self.assertTrue(matches[0][1][1] == correct[0][1][1])
        self.assertTrue(matches[1][1][1] == correct[1][1][1])
        self.assertTrue(matches[2][1][1] == correct[2][1][1])
        self.assertTrue(matches[3][1][1] == correct[3][1][1])

unittest.main()
