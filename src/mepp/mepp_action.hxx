///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_MEPP_ACTION
#define HEADER_MEPP_ACTION

#include <QObject>

#include "./Polyhedron/polyhedron.h"

class mepp_action : public QObject
{
	Q_OBJECT

	public:
		/*mepp_action( QObject *parent=0, char *name=0 ) : QObject( parent, name )
		{
		}*/

		void doSendParams(QString title, QString typeFiles, typeFuncOpenSave f)
		{
			emit(sendParams(title, typeFiles, f));
		}

		void doSendParamsOpen(QString title, QString typeFiles, int loadType, typeFuncOpenSave f)
		{
			emit(sendParamsOpen(title, typeFiles, loadType, f));
		}

	signals:
		void sendParams(QString, QString, typeFuncOpenSave);
		void sendParamsOpen(QString, QString, int, typeFuncOpenSave);
};

#endif