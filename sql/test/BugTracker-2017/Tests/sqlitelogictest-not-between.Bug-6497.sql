CREATE TABLE tab0(col0 INTEGER, col1 INTEGER, col2 INTEGER);
INSERT INTO tab0 VALUES(97,1,99);
INSERT INTO tab0 VALUES(15,81,47);
INSERT INTO tab0 VALUES(87,21,10);
SELECT * FROM tab0 WHERE + col1 NOT BETWEEN col2 + - 99 AND + 20;
DROP TABLE tab0;