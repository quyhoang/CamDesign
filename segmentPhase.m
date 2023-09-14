function segmentPhase = segmentPhase(P1,P2)
% return the angle of the vector from point P1 to point P2

complexNum = P2(1) - P1(1) + (P2(2) - P1(2))*1i;
segmentPhase = rad2deg(angle(complexNum));
end