

#make sure that Fiji.app\jars contains sqlite-jdbc-xxxx.jar
from com.ziclix.python.sql import zxJDBC


def getconn():
	db_path = r"C:\port\hudrive\OneDrive - hu-berlin.de\lab\sfv\database\data\_sfv_.db"
	DB_URL = "jdbc:sqlite:" + db_path
	DRIVER = "org.sqlite.JDBC"
	return zxJDBC.connect(DB_URL, "", "", DRIVER)


def getdata():
	conn = getconn()
	cursor = conn.cursor()
	cursor.execute("SELECT name,type FROM sqlite_master order by type;")
	data = [{"name": row[0], "type": row[1]} for row in cursor.fetchall()]
	cursor.close()
	conn.close()
	return data


data = getdata()
print data[32]
