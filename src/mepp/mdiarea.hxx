/*!
 * \file mdiarea.hxx
 * \brief MdiArea file.
 * \author Martial TOLA, CNRS-Lyon, LIRIS UMR 5205
 * \date 2010
 */ 
#ifndef HEADER_MDIAREA
#define HEADER_MDIAREA

#ifndef _MSC_VER
#pragma GCC diagnostic ignored "-Wuninitialized"
#endif
#include <QMdiArea>
#ifndef _MSC_VER
#pragma GCC diagnostic warning "-Wuninitialized"
#endif

#include <QDragEnterEvent>

class mainwindow;

enum { bNone, bLeft, bRight };

/*! 
 * \class MdiArea
 * \brief MdiArea class.
 */
class MdiArea : public QMdiArea
{
	Q_OBJECT

	public:
		/*!
		 * \brief Constructor.
		 *
		 * \param parent : parent window.
		 */
		MdiArea(QWidget *parent);
		/*!
		 * \fn setMainWindow(mainwindow* mw)
		 * \brief Set mainwindow pointer to mdiarea.
		 *
		 * \param mw mainwindow pointer.
		 */
		void setMainWindow(mainwindow* mw) { m_mw = mw; }

	protected:
		/*!
		 * \fn dragEnterEvent(QDragEnterEvent *event)
		 * \brief Accept a drop mesh file in this zone.
		 *
		 * \param event an event.
		 */
		void dragEnterEvent(QDragEnterEvent *event);
		/*!
		 * \fn dropEvent(QDropEvent *event)
		 * \brief Load drop mesh file.
		 *
		 * \param event an event.
		 */
		void dropEvent(QDropEvent *event);

		/*!
		 * \fn paintEvent(QPaintEvent *)
		 * \brief paint Mepp mdiarea.
		 *
		 * \param paintEvent an event.
		 */
		void paintEvent(QPaintEvent *);

	private:
		mainwindow* m_mw;	//!< mainwindow
		int bType;			//!< mouse bouton type (bNone, bLeft, bRight)
};

#endif
