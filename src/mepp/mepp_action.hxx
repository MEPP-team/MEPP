/*!
 * \file mepp_action.hxx
 * \brief mepp_action file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */
#ifndef HEADER_MEPP_ACTION
#define HEADER_MEPP_ACTION

#include <QObject>

#include "./Polyhedron/polyhedron.h"

/*! 
 * \class mepp_action
 * \brief mepp_action class.
 */
class mepp_action : public QObject
{
	Q_OBJECT

	public:
		/*mepp_action( QObject *parent=0, char *name=0 ) : QObject( parent, name )
		{
		}*/

		/*!
		 * \fn doSendParams(QString title, QString typeFiles, typeFuncOpenSave f)
		 * \brief Emit a signal for loading or saving a mesh from a specific component.
		 *
		 * \param title title of the dialog box, for example: tr("Open Mesh File(s) - from CGAL_Example")
		 * \param typeFiles title of the dialog box, for example: tr("OFF files (*.off)")
 		 * \param f pointer to a typeFuncOpenSave function.
		 */
		void doSendParams(QString title, QString typeFiles, typeFuncOpenSave f)
		{
			emit(sendParams(title, typeFiles, f));
		}

		/*!
		 * \fn doSendParamsOpen(QString title, QString typeFiles, int loadType, typeFuncOpenSave f)
		 * \brief Emit a signal for loading or saving a mesh from a specific component.
		 *
		 * \param title title of the dialog box, for example: tr("Open Mesh File(s) - from CGAL_Example")
		 * \param typeFiles title of the dialog box, for example: tr("OFF files (*.off)")
 		 * \param loadType (Normal, Space or Time).
 		 * \param f pointer to a typeFuncOpenSave function.
		 */
		void doSendParamsOpen(QString title, QString typeFiles, int loadType, typeFuncOpenSave f)
		{
			emit(sendParamsOpen(title, typeFiles, loadType, f));
		}

	signals:
		void sendParams(QString, QString, typeFuncOpenSave);
		void sendParamsOpen(QString, QString, int, typeFuncOpenSave);
};

#endif