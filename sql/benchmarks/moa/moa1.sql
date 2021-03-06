 -- show tables like '_t%';
 SELECT * FROM MOA_DD WHERE name='Video' AND type='match';
 SELECT * FROM MOA_DD WHERE name='Creation' AND type='match';
 SELECT * FROM MOA_DD WHERE name='CreationCreator' AND type='match';
 SELECT * FROM MOA_DD WHERE name='CreationCreator' AND type='primary';
 SELECT * FROM MOA_DD WHERE name='Creator' AND type='match';
 SELECT * FROM MOA_DD WHERE name='Individual' AND type='match';
 create table _t1 (CreationPK int, Title varchar(255), CreationDate varchar(255));
 INSERT INTO _t1 SELECT * from Creation WHERE (Creation.Title = 'Titanic');
 create table _t2 (ID int, IMDB varchar(255), Image varchar(255), MediaInformationFK int);
 INSERT INTO _t2 SELECT Video.ID, Video.IMDB, Video.Image, Video.MediaInformationFK FROM Video,_t1 WHERE (Video.ID = _t1.CreationPK);
 create table _t3 (CreationFK int, CreatorFK int);
 INSERT INTO _t3 SELECT CreationCreator.CreationFK, CreationCreator.CreatorFK FROM CreationCreator,_t1 WHERE (CreationCreator.CreationFK = _t1.CreationPK);
 create table _t4 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t4 SELECT Creator.CreatorPK, Creator.Role, Creator.IndividualFK FROM Creator,_t3 WHERE (Creator.CreatorPK = _t3.CreatorFK);
 create table _t5 (IndividualPK int, PersonName varchar(255));
 INSERT INTO _t5 SELECT Individual.IndividualPK, Individual.PersonName FROM Individual,_t4 WHERE (Individual.IndividualPK = _t4.IndividualFK);
 SELECT * from _t2;
 SELECT * from _t1;
 SELECT * from _t3;
 SELECT * from _t4;
 SELECT * from _t5;
 SELECT * from Video;
 SELECT * from Creation;
 SELECT * from CreationCreator;
 SELECT * from Creator;
 SELECT * from Individual;
 create table _t6 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t6 SELECT * from Creator WHERE (Creator.Role = 'Director');
 create table _t7 (CreationFK int, CreatorFK int);
 INSERT INTO _t7 SELECT CreationCreator.CreationFK, CreationCreator.CreatorFK FROM CreationCreator,_t6 WHERE (CreationCreator.CreatorFK = _t6.CreatorPK);
 create table _t8 (IndividualPK int, PersonName varchar(255));
 INSERT INTO _t8 SELECT Individual.IndividualPK, Individual.PersonName FROM Individual,_t6 WHERE (Individual.IndividualPK = _t6.IndividualFK);
 SELECT * from _t7;
 SELECT * from _t6;
 SELECT * from _t8;
 create table _t9 (CreationPK int, Title varchar(255), CreationDate varchar(255));
 INSERT INTO _t9 SELECT * from Creation WHERE (Creation.CreationDate = '1996');
 create table _t10 (ID int, IMDB varchar(255), Image varchar(255), MediaInformationFK int);
 INSERT INTO _t10 SELECT Video.ID, Video.IMDB, Video.Image, Video.MediaInformationFK FROM Video,_t9 WHERE (Video.ID = _t9.CreationPK);
 create table _t11 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t11 SELECT * from Creator WHERE (Creator.Role = 'Writer');
 create table _t12 (CreationFK int, CreatorFK int);
 INSERT INTO _t12 SELECT CreationCreator.CreationFK, CreationCreator.CreatorFK FROM CreationCreator,_t11 WHERE (CreationCreator.CreatorFK = _t11.CreatorPK);
 create table _t13 (CreationFK int, CreatorFK int);
 INSERT INTO _t13 SELECT _t12.CreationFK, _t12.CreatorFK FROM _t12,_t9 WHERE (_t12.CreationFK = _t9.CreationPK);
 create table _t14 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t14 SELECT _t11.CreatorPK, _t11.Role, _t11.IndividualFK FROM _t11,_t13 WHERE (_t11.CreatorPK = _t13.CreatorFK);
 create table _t15 (IndividualPK int, PersonName varchar(255));
 INSERT INTO _t15 SELECT Individual.IndividualPK, Individual.PersonName FROM Individual,_t11 WHERE (Individual.IndividualPK = _t11.IndividualFK);
 create table _t16 (IndividualPK int, PersonName varchar(255));
 INSERT INTO _t16 SELECT _t15.IndividualPK, _t15.PersonName FROM _t15,_t14 WHERE (_t15.IndividualPK = _t14.IndividualFK);
 SELECT * from _t10;
 SELECT * from _t9;
 SELECT * from _t13;
 SELECT * from _t14;
 SELECT * from _t16;
 create table _t17 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t17 SELECT * from Creator WHERE (Creator.Role = 'Director');
 create table _t18 (CreationFK int, CreatorFK int);
 INSERT INTO _t18 SELECT CreationCreator.CreationFK, CreationCreator.CreatorFK FROM CreationCreator,_t17 WHERE (CreationCreator.CreatorFK = _t17.CreatorPK);
 create table _t19 (CreationPK int, Title varchar(255), CreationDate varchar(255));
 INSERT INTO _t19 SELECT Creation.CreationPK, Creation.Title, Creation.CreationDate FROM Creation,_t18 WHERE (Creation.CreationPK = _t18.CreationFK);
 SELECT * from _t17;
 SELECT * from _t18;
 SELECT * from _t19;
 create table _t20 (IndividualPK int, PersonName varchar(255));
 INSERT INTO _t20 SELECT * from Individual WHERE (Individual.PersonName = 'Spielberg,+Steven');
 create table _t21 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t21 SELECT * from Creator WHERE (Creator.Role = 'Director');
 create table _t22 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t22 SELECT _t21.CreatorPK, _t21.Role, _t21.IndividualFK FROM _t21,_t20 WHERE (_t21.IndividualFK = _t20.IndividualPK);
 create table _t23 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t23 SELECT _t21.CreatorPK, _t21.Role, _t21.IndividualFK FROM _t21,_t22 WHERE (_t21.CreatorPK = _t22.CreatorPK);
 create table _t24 (CreationFK int, CreatorFK int);
 INSERT INTO _t24 SELECT CreationCreator.CreationFK, CreationCreator.CreatorFK FROM CreationCreator,_t21 WHERE (CreationCreator.CreatorFK = _t21.CreatorPK);
 create table _t25 (CreationFK int, CreatorFK int);
 INSERT INTO _t25 SELECT _t24.CreationFK, _t24.CreatorFK FROM _t24,_t22 WHERE (_t24.CreatorFK = _t22.CreatorPK);
 create table _t26 (CreationPK int, Title varchar(255), CreationDate varchar(255));
 INSERT INTO _t26 SELECT Creation.CreationPK, Creation.Title, Creation.CreationDate FROM Creation,_t24 WHERE (Creation.CreationPK = _t24.CreationFK);
 create table _t27 (CreationPK int, Title varchar(255), CreationDate varchar(255));
 INSERT INTO _t27 SELECT _t26.CreationPK, _t26.Title, _t26.CreationDate FROM _t26,_t25 WHERE (_t26.CreationPK = _t25.CreationFK);
 SELECT * from _t20;
 SELECT * from _t22;
 SELECT * from _t23;
 SELECT * from _t25;
 SELECT * from _t27;
 create table _t29 (IndividualPK int, PersonName varchar(255));
 INSERT INTO _t29 SELECT * from Individual WHERE (Individual.PersonName LIKE '%Spielberg%');
 create table _t28 (CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t28 SELECT * from Creator WHERE (Creator.Role = 'Director');
 create table _t30 (IndividualPK int, PersonName varchar(255), CreatorPK int, Role varchar(255), IndividualFK int);
 INSERT INTO _t30 SELECT _t29.IndividualPK, _t29.PersonName, _t28.CreatorPK, _t28.Role, _t28.IndividualFK FROM _t29,_t28 WHERE (_t29.IndividualPK = _t28.IndividualFK);
 create table _t31 (IndividualPK int, PersonName varchar(255), CreatorPK int, Role varchar(255), IndividualFK int, CreationFK int, CreatorFK int);
 INSERT INTO _t31 SELECT _t30.IndividualPK, _t30.PersonName, _t30.CreatorPK, _t30.Role, _t30.IndividualFK, CreationCreator.CreationFK, CreationCreator.CreatorFK FROM _t30,CreationCreator WHERE (_t30.CreatorPK = CreationCreator.CreatorFK);
 create table _t32 (IndividualPK int, PersonName varchar(255), CreatorPK int, Role varchar(255), IndividualFK int, CreationFK int, CreatorFK int, CreationPK int, Title varchar(255), CreationDate varchar(255));
 INSERT INTO _t32 SELECT _t31.IndividualPK, _t31.PersonName, _t31.CreatorPK, _t31.Role, _t31.IndividualFK, _t31.CreationFK, _t31.CreatorFK, Creation.CreationPK, Creation.Title, Creation.CreationDate FROM _t31,Creation WHERE (_t31.CreationFK = Creation.CreationPK);
 SELECT * from _t32;
 drop table _t1;
 drop table _t2;
 drop table _t3;
 drop table _t4;
 drop table _t5;
 drop table _t6;
 drop table _t7;
 drop table _t8;
 drop table _t9;
 drop table _t10;
 drop table _t11;
 drop table _t12;
 drop table _t13;
 drop table _t14;
 drop table _t15;
 drop table _t16;
 drop table _t17;
 drop table _t18;
 drop table _t19;
 drop table _t20;
 drop table _t21;
 drop table _t22;
 drop table _t23;
 drop table _t24;
 drop table _t25;
 drop table _t26;
 drop table _t27;
 drop table _t28;
 drop table _t29;
 drop table _t30;
 drop table _t31;
 drop table _t32;
