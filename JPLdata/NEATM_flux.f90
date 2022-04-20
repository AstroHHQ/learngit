
program main
implicit none
!--------------------------------------
!constants:
integer(4),parameter :: epn = 2  !number of epochs
real(8),parameter :: epsi = 0.9
real(8),parameter :: sigmas = 5.67d-8
real(8),parameter :: f_solar = 1367.5
real(8),parameter :: h = 6.62607015d-34
real(8),parameter :: c = 3d8 
real(8),parameter :: kb = 1.380649d-23
real(8),parameter :: au = 1.496d11
real(8),parameter :: pi = 3.141592653589793
real(8),parameter :: SR = 6.9550826d8
!--------------------------------------
!variables:
character (len=15) ::syita,sDia,sdao,sdast,salpha,swlenth,sA
real(8) :: G,pv,qph,A,yita,Dia,wlenth,flux_i,alpha,dast,dao
!slope parameter,geometric albedo, phase integral, bond albedo, beaming parameter, diameter
!real(8) :: dast,dao,dobs,alpha(epn),cos_alpha,astp(3),obsp(3)
! heliocentric distance, ast-observer distance, observer-sun ditance, phase angle, cosine of phase angle, ast cartesian coord, obs coord
integer(4) :: i,j,k
real(8),external :: dist
call get_command_argument(1,sdast)
read(sdast,"(f8.5)")dast
call get_command_argument(2,sdao)
read(sdao,"(f8.5)")dao
call get_command_argument(3,salpha)
read(salpha,"(f8.5)")alpha
call get_command_argument(4,sDia)
read(sDia,"(f8.5)")Dia
call get_command_argument(5,swlenth)
read(swlenth,"(f8.5)")wlenth
call get_command_argument(6,syita)
read(syita,"(f8.5)")yita
call get_command_argument(7,sA)
read(sA,"(f8.5)")A

!print
print *,' hello' 
print 100
100 format (5x,'dast:', 5x, 'dao:', 5x, 'alpha:', 5x, 'dia:',5x, 'wlenth:',5x, 'yita:',5x, 'A:')
print 200, dast,dao,alpha,dia,wlenth,yita,A 
200 format(f10.6,f10.6,f10.6,f10.6,f10.6,f10.6,f10.6)

!call get_command_argument(8,the_name)

call neatm_flux(flux_i,dast,dao,alpha,Dia,wlenth,yita,A,h,c,epsi,f_solar,sigmas,au,kb)
print 300
300 format (7x,'h:', 7x, 'c:', 7x, 'epsi:', 7x, 'F_un:',7x, 'sigmas:',7x, 'au:',7x, 'kb')
print 400, h,c,epsi,f_solar,sigmas,au,kb 
400 format(e10.3,e10.3,e10.3,e10.3,e10.3,e10.3,e10.3)

print'(/,"Flux_i = ",e10.4)',flux_i

end program



!--------------------------------------
!calculate the theoretical thermal flux in mjy
subroutine neatm_flux(flux,dast,dao,alpha,Dia,wlenth,yita,A,h,c,epsi,f_solar,sigmas,au,kb)
implicit none
real(8),parameter :: pi = 3.141592653589793
integer(4),parameter :: dp = 100
real(8) :: T_ss,phi,theta,dphi,dtheta,temp(dp,dp)
real(8) :: dast,dao,dobs,alpha,flux,Dia,wlenth,yita,A,h,c,epsi,f_solar,sigmas,au,flux_con,kb
integer(4) :: i,j,k,nj

dphi = pi/dp
dtheta = pi/dp
T_ss = ((1 - A) * f_solar / epsi / yita / sigmas / dast ** 2) ** 0.25 !subsolar temperature
nj = int(floor(((pi/2 - alpha) + pi / 2) / (pi / dp)))
temp = 1d0
flux_con = epsi * Dia ** 2 * pi  * h * c ** 2  / (wlenth ** 5)

do j = 1,dp
	do k = nj,dp
		theta = -pi / 2 + (k - 1) * dtheta
		phi = -pi / 2 + (j - 1) * dphi
		temp(j,k) = T_ss * abs(cos(theta)) ** 0.25 * abs(cos(phi)) ** 0.25
		flux = flux + flux_con / 4 * abs(cos(alpha)) * abs(cos(alpha - theta)) / (2 * (dao * au) ** 2) &
		/ (exp(h * c / (wlenth * kb * temp(j,k))) - 1) * dphi * dtheta * wlenth ** 2 / c * 1d29 ! obtain flux in unit of mjy
	end do
end do



end subroutine



!--------------------------------------
!calculate the distance between two points
function dist(x1,y1,z1,x2,y2,z2)
implicit none
real(8) dist,x1,y1,z1,x2,y2,z2

dist = sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

end function







