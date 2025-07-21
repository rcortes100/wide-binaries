-- Basado en el query de ElBadry 2018, pero optimizado
-- 15 de abril de 2025
-- Para búsqueda automatica, por ejemplo en Python, hay que cambiar el intervalo de
-- búsqueda dentro del healpix en las tablas pri y sec
-- Declaracion de tablas pri y sec (deben ser identicas)
-- Notas sobre Healpix:
-- Healpix_index van de 0 a N-1, donde para cada nivel:
-- N = 12*4^nivel
-- nivel 0: 12
-- nivel 1: 48
-- nivel 2: 192
-- nivel 3: 768
-- etc.

WITH pri AS (
  SELECT *
  FROM gaiadr2.gaia_source
  WHERE parallax BETWEEN 5 AND 1000  	   	-- Estrellas entre 1 y 200 pc
    AND GAIA_HEALPIX_INDEX(2,source_id)
    	BETWEEN 185 AND 185			-- NÃºmero de healpix
    AND parallax_over_error > 20     	   	-- Error en paralaje menor al 5%
    AND phot_g_mean_flux_over_error > 50   	-- Error en brillo G menor al 2%
    AND phot_rp_mean_flux_over_error > 20  	-- Error en brillo RP menor al 5%
    AND phot_bp_mean_flux_over_error > 20  	-- Error en brillo BP menor al 5%
	),
  sec AS (
  SELECT *
  FROM gaiadr2.gaia_source
  WHERE parallax BETWEEN 5 AND 1000		-- Estrellas entre 1 y 200 pc
    AND GAIA_HEALPIX_INDEX(2,source_id)
    	BETWEEN 185 AND 185			-- NÃºmero de healpix
    AND parallax_over_error > 20    		-- Error en paralaje menor al 5%
    AND phot_g_mean_flux_over_error > 50   	-- Error en brillo G menor al 2%
    AND phot_rp_mean_flux_over_error > 20	-- Error en brillo RP menor al 5%
    AND phot_bp_mean_flux_over_error > 20	-- Error en brillo BP menor al 5%
	)

-- Seleccion de tabla final
SELECT
    pri.source_id AS source_id_1, sec.source_id AS source_id_2,					-- Source identifier
    pri.ra AS ra_1, pri.dec AS dec_1, sec.ra AS ra_2, sec.dec AS dec_2, 	               	-- Coordinates
--	GAIA_HEALPIX_INDEX(2,pri.source_id) AS hp_1,
--	CAST (pri.source_id/576460752303423488 AS INTEGER) AS my_hp_1,
--	GAIA_HEALPIX_INDEX(2,sec.source_id) AS hp_2,
--	CAST (sec.source_id/576460752303423488 AS INTEGER) AS my_hp_2,
    pri.ra_error AS ra_err_1, pri.dec_error AS dec_err_1,					-- Coordinate error primary
    sec.ra_error AS ra_err_2, sec.dec_error AS dec_err_2,					-- Coordinate error secondary
    pri.parallax AS parallax_1, sec.parallax AS parallax_2,					-- Parallax
    pri.parallax_over_error AS par_over_err_1, sec.parallax_over_error AS par_over_err_2,	-- Parallax error
    pri.pmra AS pmra_1, pri.pmdec AS pmdec_1, sec.pmra AS pmra_2, sec.pmdec AS pmdec_2,		-- Proper motion
    pri.pmra_error AS pmra_err_1, pri.pmdec_error AS pmdec_err_1,				-- Proper motion error primary
    sec.pmra_error AS pmra_err_2, sec.pmdec_error AS pmdec_err_2,				-- Proper motion error secondary
    pri.radial_velocity AS rad_vel_1, sec.radial_velocity AS rad_vel_2,				-- Radial velocity
    pri.radial_velocity_error AS rad_vel_err_1,							-- Radial velocity error pri
    sec.radial_velocity_error AS rad_vel_err_2,							-- Radial velocity error sec
    pri.phot_g_mean_mag AS phot_g_1, sec.phot_g_mean_mag AS phot_g_2,				-- G mag
    pri.phot_g_mean_flux_over_error AS phot_g_over_err_1,					-- G mag error primary
    sec.phot_g_mean_flux_over_error AS phot_g_over_err_2,					-- G mag error secondary
    pri.phot_bp_mean_mag AS phot_bp_1, sec.phot_bp_mean_mag AS phot_bp_2,			-- BP mag
    pri.phot_bp_mean_flux_over_error AS phot_bp_over_err_1,					-- BP mag error primary
    sec.phot_bp_mean_flux_over_error AS phot_bp_over_err_2,					-- BP mag error secondary
    pri.phot_rp_mean_mag AS phot_rp_1, sec.phot_rp_mean_mag AS phot_rp_2,			-- RP mag
    pri.phot_rp_mean_flux_over_error AS phot_rp_1,						-- RP mag error primary
    sec.phot_rp_mean_flux_over_error AS phot_rp_2,						-- RP mag error secondary
    pri.astrometric_chi2_al AS astro_chi2_al_1,							-- Astrometric chi2 factor pri
    sec.astrometric_chi2_al AS astro_chi2_al_2,							-- Astrometric chi2 factor sec
    pri.astrometric_chi2_al AS astro_chi2_al_1,							-- Astrometric chi2 factor pri
    sec.astrometric_chi2_al AS astro_chi2_al_2,							-- Astrometric chi2 factor sec
    pri.astrometric_n_good_obs_al AS astro_n_goodobs_al_1,					-- Astrometric good obs pri
    sec.astrometric_n_good_obs_al AS astro_n_goodobs_al_2,					-- Astrometric good obs sec
    pri.phot_bp_rp_excess_factor AS phot_bp_rp_xs_1,						-- Fotometric excess factor pri
    sec.phot_bp_rp_excess_factor AS phot_bp_rp_xs_2,						-- Fotometric excess factor sec
    pri.rv_nb_transits AS rv_nb_tran_1,								-- RV transits pri
    sec.rv_nb_transits AS rv_nb_tran_2,								-- RV transits sec
    DISTANCE(POINT('ICRS', pri.ra, pri.dec),
             POINT('ICRS', sec.ra, sec.dec)) AS pairdistance,					-- Pair distance
    SQRT(POWER(pri.ra - sec.ra, 2) + POWER(pri.dec - sec.dec, 2)) AS angular_distance		-- Angular Distance

-- Join de las tablas	
FROM pri
JOIN sec ON
    pri.source_id < sec.source_id  			                             		-- Estrellas distintas
--    pri.phot_g_mean_mag > sec.phot_g_mean_mag							-- Ordenamiento por brillo
    AND 1 = CONTAINS(POINT('ICRS', sec.ra, sec.dec),						-- Distancia angular dentro de un radio
                 CIRCLE('ICRS', pri.ra, pri.dec,						--  de 50000AU: theta/arcsec <= 50 * parallax/mas
                 1.4e-2*pri.parallax))								-- 
    AND ABS(1/pri.parallax - 1/sec.parallax) -							-- Distancia entre las estrellas debe ser
          2*0.01745*DISTANCE(POINT('ICRS', pri.ra, pri.dec),					-- menor al doble de su separaciÃ³n dentro de tres
          		     POINT('ICRS', sec.ra, sec.dec))/					-- desviaciones estÃ¡ndar:
            pri.parallax <									-- Delta d - 2s <= 3 sigma delta d
        3*SQRT(POWER(pri.parallax_error,2)/POWER(pri.parallax, 4) +				--
               POWER(sec.parallax_error,2)/POWER(sec.parallax, 4))				--
    AND SQRT(POWER((pri.pmra - sec.pmra), 2) +							-- Diferencia entre movimientos propios
    	     POWER((pri.pmdec - sec.pmdec), 2))							-- debe ser menor a la diferencia de mov. prop.
     	- (7.42e-3 * POWER(pri.parallax, 1.5)							-- orbital  mas tres veces la desviacion estandar:
     	   *POWER(DISTANCE(POINT('ICRS', pri.ra, pri.dec),					-- Delta mu <= Delta mu_orb + 3 sigma Delta mu
			   POINT('ICRS', sec.ra, sec.dec)), -0.5)) <				--
	3*sqrt(((POWER((pri.pmra_error), 2)							--
		 +POWER((sec.pmra_error), 2)) * POWER((pri.pmra - sec.pmra), 2)			--
		 +(POWER((pri.pmdec_error), 2) + POWER((sec.pmdec_error), 2)) *			--
		  POWER((pri.pmdec - sec.pmdec), 2))/						--
		(POWER((pri.pmra - sec.pmra), 2) + POWER((pri.pmdec - sec.pmdec), 2))		--
		)										-- 
    AND SQRT(((POWER((pri.pmra_error), 2) +							-- Movimiento propio entre primaria y secundaria
	     POWER((sec.pmra_error), 2)) * POWER(						-- dentro de 1.5 veces el error ?
	     (pri.pmra - sec.pmra), 2) + (POWER(						--
	     (pri.pmdec_error), 2) + POWER(							--
	     (sec.pmdec_error), 2)) * POWER((pri.pmdec -					--
	     sec.pmdec), 2))/(POWER((pri.pmra - sec.pmra),					--
	     2) + POWER((pri.pmdec - sec.pmdec), 2))) < 1.5					--
ORDER BY source_id_1 ASC;
 

 --   AND ABS(pri.pmra - sec.pmra) < 1  -- Diferencia en movimiento propio en ascension recta
 --   AND ABS(pri.pmdec - sec.pmdec) < 1  -- Diferencia en movimiento propio en declinacion
 --    AND SQRT(POWER(pri.ra - sec.ra, 2) + POWER(pri.dec - sec.dec, 2)) < 0.1  -- Distancia angular en grados
 --   AND ABS(pri.parallax - sec.parallax) < 0.5 ; -- Diferencia en paralaje para asegurar proximidad real
