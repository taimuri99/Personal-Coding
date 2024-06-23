# Hamiltonian: Harmonic potential
#
def hamiltonian(r):
    pot = 0.5 * k * (r - r_eq)**2 
    return pot
#
# Monte Carlo code	
#
acc = 0
mc_dst = 0.0
mc_dst2 = 0.0
mc_ener = 0.0

r = r_ini

for point in range(1,points):
   r_new = r + step * ( random() - 0.5 ) * 2.0
   v     = hamiltonian(r)
   v_new = hamiltonian(r_new)
   if v_new < v:				# Downhill move always accepted
      r = r_new
      v = v_new
      acc = acc + 1
   else:					# Uphill move accepted with luck
   
      A_mov = exp( -beta *(v_new - v) )
      if A_mov > random():
         r = r_new
         v = v_new
         acc = acc + 1
						# Update regardless of acceptance!
   mc_dst  = mc_dst + r
   mc_dst2 = mc_dst2 + r*r
   mc_ener = mc_ener + v

#
# Print out averages
#
mc_dst_av  = mc_dst / float(points)
mc_dst2_av = mc_dst2 / float(points)
mc_ener_av = mc_ener / float(points)

print 'Acceptance is:', 100 * acc / float(points), '%'
print ' '
print 'Monte Carlo Averages:'
print 'Mean distance, Ang:                 ', mc_dst_av
print 'Mean squared distance, Ang**2       ', mc_dst2_av
print 'Mean MC square displacemet, Ang**2  ', mc_dst2_av - mc_dst_av*mc_dst_av
print 'Mean MC potential energy, kcal/mol: ', mc_ener_av
potential energy, kcal/mol: ', mc_ener_av
  	
