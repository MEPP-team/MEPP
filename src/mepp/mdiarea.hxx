///////////////////////////////////////////////////////////////////////////
// Author: Martial TOLA
// Year: 2010
// CNRS-Lyon, LIRIS UMR 5205
/////////////////////////////////////////////////////////////////////////// 
#ifndef HEADER_MDIAREA
#define HEADER_MDIAREA

#include <QMdiArea>
#include <QDragEnterEvent>

class mainwindow;

class MdiArea : public QMdiArea
{
	Q_OBJECT

	public:
		MdiArea(QWidget *parent);
		void setMainWindow(mainwindow* mw) { m_mw = mw; }

	protected:
		void dragEnterEvent(QDragEnterEvent *event);
		void dropEvent(QDropEvent *event);

	private:
		mainwindow* m_mw;
};

#endif
