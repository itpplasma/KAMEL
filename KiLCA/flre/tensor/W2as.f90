!complex(dpc), dimension(0:1, 0:3, -Nmax:Nmax) :: W2

W2(0,0,l) =  &
((-I)*sqrt2p*nu(1)*Vt(1)*(-(dom(2)*nu(2)) + 2*nu(4) + 5*kp(2)*nu(2)*Vt(2) -  &
   (3*I)*dom(1)*(nu(3) + kp(2)*nu(1)*Vt(2)) + (2 + F11m)*kp(4)*Vt(4)))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(0,1,l) =  &
((-I)*sqrt2p*kp(1)*Vt(3)*(-(dom(1)*nu(1)*(nu(2) - (-1 + F11m)*kp(2)*Vt(2))) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(0,2,l) =  &
(I*sqrt2p*(dom(1) + I*nu(1))*Vt(3)* &
  (dom(1)*nu(1)*(nu(2) - (-1 + F11m)*kp(2)*Vt(2)) +  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(0,3,l) =  &
(sqrt2p*kp(1)*((-1 + F11m)*dom(3)*nu(1) - (1 + 2*F11m)*dom(1)*nu(3) +  &
   I*dom(2)*(3*(-1 + F11m)*nu(2) - kp(2)*Vt(2)) -  &
   I*(6*nu(4) + (5 + 2*F11m)*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4)))*Vt(5))/ &
 (I*dom(3)*nu(3) + 2*kp(2)*nu(4)*Vt(2) -  &
  3*dom(2)*(nu(4) + kp(2)*nu(2)*Vt(2)) + 3*kp(4)*nu(2)*Vt(4) -  &
  I*dom(1)*(2*nu(5) + 6*kp(2)*nu(3)*Vt(2) + 3*kp(4)*nu(1)*Vt(4)) +  &
  kp(6)*Vt(6))

W2(1,0,l) =  &
((-I)*sqrt2p*kp(1)*Vt(3)*(-(dom(1)*nu(1)*(nu(2) - (-1 + F11m)*kp(2)*Vt(2))) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(1,1,l) =  &
((-I)*sqrt2p*dom(1)*Vt(3)* &
  (-(dom(1)*nu(1)*(nu(2) - (-1 + F11m)*kp(2)*Vt(2))) -  &
   I*(2*nu(4) + 3*kp(2)*nu(2)*Vt(2) + kp(4)*Vt(4))))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))

W2(1,2,l) =  &
-((sqrt2p*kp(1)*(I*(-1 + F11m)*dom(3)*nu(1) + 2*nu(4) + 3*kp(2)*nu(2)*Vt(2) +  &
    dom(2)*(-((-1 + F11m)*nu(2)) + kp(2)*Vt(2)) -  &
    I*dom(1)*(3*nu(3) + 2*kp(2)*nu(1)*Vt(2)) + kp(4)*Vt(4))*Vt(5))/ &
  ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
   (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2)))))

W2(1,3,l) =  &
(sqrt2p*dom(1)*((-I)*(-1 + F11m)*dom(3)*nu(1) + I*(1 + 2*F11m)*dom(1)*nu(3) -  &
   6*nu(4) - (5 + 2*F11m)*kp(2)*nu(2)*Vt(2) +  &
   dom(2)*(3*(-1 + F11m)*nu(2) - kp(2)*Vt(2)) - kp(4)*Vt(4))*Vt(5))/ &
 ((dom(1)*nu(1) + I*kp(2)*Vt(2))*(dom(1)*nu(1) + I*(nu(2) + kp(2)*Vt(2)))* &
  (dom(1)*nu(1) + I*(2*nu(2) + kp(2)*Vt(2))))
