import lmfdb

pqs = []
for p in [2,3,5,7,13]:
#  print("p =", p, "\n")
  for q in prime_range(2,floor(10000/p)):
    if q == p:
      continue
    newforms = lmfdb.db.mf_newforms.search({"weight":2,"char_order":1,"level":p*q,"analytic_rank":0,"analytic_rank_proved":True,"atkin_lehner_eigenvals":sorted([[p,-1],[q,1]])})
    if len(list(newforms)) != 0:
      pass
      #print("q =", q, "is OK!")
    else:
      print("p = ", p, ", q =", q, ": have to look for oldspace")
      pqs.append([p,q])
#  print("\n")

for pq in pqs:
  p, q = pq[0], pq[1]
  for k in range(1,10):
    newforms = lmfdb.db.mf_newforms.search({"weight":2,"char_order":1,"level":p^k*q,"analytic_rank":0,"analytic_rank_proved":True,"atkin_lehner_eigenvals":sorted([[p,-1],[q,1]])})
    if len(list(newforms)) != 0:
      print("p =", p, ", q =", q, ": k_min =", k)
      break
  for l in range(1,10):
    newforms = lmfdb.db.mf_newforms.search({"weight":2,"char_order":1,"level":p*q^l,"analytic_rank":0,"analytic_rank_proved":True,"atkin_lehner_eigenvals":sorted([[p,-1],[q,1]])})
    if len(list(newforms)) != 0:
      print("p =", p, ", q =", q, ": l_min =", l)
      break

p = 5
q = 2
k = 2
newforms = lmfdb.db.mf_newforms.search({"weight":2,"char_order":1,"level":p^k*q,"analytic_rank":0,"analytic_rank_proved":True,"atkin_lehner_eigenvals":sorted([[p,-1],[q,1]])})
if len(list(newforms)) != 0:
  print("p =", p, ", q =", q, ": k_min =", k)
else:
  print("not found")