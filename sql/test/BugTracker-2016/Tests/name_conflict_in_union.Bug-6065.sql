WITH t0(so) AS (
	SELECT ROW_NUMBER() OVER (ORDER BY v0 ASC) AS so
FROM (values (0)) as t(v0))
SELECT 1 AS k1, a3.so AS k2 FROM t0 AS a3
UNION ALL
SELECT 2 AS k2, a8.so AS k2 FROM t0 AS a8;