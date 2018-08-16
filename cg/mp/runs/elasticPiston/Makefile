#
#   Makefile for the elastic piston
#


ogen = $(Overture)/bin/ogen
# ogen = /home/dan/Documents/Overture.sf3/bin/ogen

grids: elasticPiston radial rotated

elasticPiston:
	$(ogen) -noplot elasticPistonGrid -factor=1
	$(ogen) -noplot elasticPistonGrid -factor=2
	$(ogen) -noplot elasticPistonGrid -factor=4
	$(ogen) -noplot elasticPistonGrid -factor=8
	$(ogen) -noplot elasticPistonGrid -factor=16
	$(ogen) -noplot elasticPistonGrid -factor=1 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=2 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=4 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=8 -per=1
	$(ogen) -noplot elasticPistonGrid -factor=16 -per=1


# Narrow grid in x-direction
elasticPistonN:
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=1 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=2 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=4 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -factor=8 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=1 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=2 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=4 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=8 
	$(ogen) -noplot elasticPistonGrid -xb=.2 -prefix=elasticPistonGridN -per=1 -factor=16

# coarse grid in x: 
elasticPistonX:
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=1 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=2 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=4 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -factor=8 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=1 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=2 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=4 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=8 
	$(ogen) -noplot elasticPistonGrid -dsx=.2 -prefix=elasticPistonGridX -per=1 -factor=16


radial:
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=1
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=2
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=4
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=8
	$(ogen) -noplot radialElasticPistonGrid -interp=e -factor=16
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=1
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=2
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=4
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=8
	$(ogen) -noplot radialTravelingWaveGrid -interp=e -factor=16

rotated:
	$(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=0
	$(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=30
	$(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=0
	# $(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=30
	# $(ogen) -noplot elasticPistonGridRotated -factor=1 -per=1 -thetad=45
	# $(ogen) -noplot elasticPistonGridRotated -factor=2 -per=1 -thetad=45
	# $(ogen) -noplot elasticPistonGridRotated -factor=4 -per=1 -thetad=45
	# $(ogen) -noplot elasticPistonGridRotated -factor=8 -per=1 -thetad=45




mp:
	cd ../../; \
	make

