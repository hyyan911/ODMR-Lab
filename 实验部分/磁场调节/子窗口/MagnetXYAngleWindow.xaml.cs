﻿using CodeHelper;
using CodeHelper.ResizeElements;
using Controls.Charts;
using MathNet.Numerics;
using ODMR_Lab.Windows;
using ODMR_Lab.实验部分.磁场调节;
using System;
using System.Collections.Generic;
using System.Data;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using Window = System.Windows.Window;

using ODMR_Lab.基本控件;
using Controls;

namespace ODMR_Lab.磁场调节
{
    /// <summary>
    /// MagnetXYZWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MagnetXYAngleWindow : Window
    {
        public ItemsList<CWPointObject> CWPoints { get; private set; } = new ItemsList<CWPointObject>();

        private NumricChartData1D DataBp = new NumricChartData1D("沿轴磁场(Gauss)", "磁场扫描数据", ChartDataType.Y);
        private NumricChartData1D DataBv = new NumricChartData1D("垂直轴磁场(Gauss)", "磁场扫描数据", ChartDataType.Y);
        private NumricChartData1D DataB = new NumricChartData1D("总磁场(Gauss)", "磁场扫描数据", ChartDataType.Y);
        private NumricChartData1D DataLoc = new NumricChartData1D("位置", "磁场扫描数据", ChartDataType.X);

        private NumricChartData1D DataCW1 = new NumricChartData1D("CW数据1", "磁场扫描CW数据", ChartDataType.Y);
        private NumricChartData1D DataCW2 = new NumricChartData1D("CW数据2", "磁场扫描CW数据", ChartDataType.Y);
        private NumricChartData1D Freq1 = new NumricChartData1D("频率1(MHz)", "磁场扫描CW数据", ChartDataType.X);
        private NumricChartData1D Freq2 = new NumricChartData1D("频率2(MHz)", "磁场扫描CW数据", ChartDataType.X);

        private CWPointObject CurrentPoint = null;

        private bool IsAngleWindow = false;

        private DisplayPage ParentPage { get; set; } = null;

        public MagnetXYAngleWindow(string windowtitle, bool isAngleWindow, DisplayPage parent)
        {
            InitializeComponent();
            ParentPage = parent;
            IsAngleWindow = isAngleWindow;
            Title = windowtitle;
            WindowTitle.Content = windowtitle;
            CodeHelper.WindowResizeHelper helper = new CodeHelper.WindowResizeHelper();
            helper.RegisterWindow(this, MinBtn, null, null, dragHeight: 30);

            Chart.DataSource.AddRange(new List<ChartData1D>()
            {
                DataBv,
                DataBp,
                DataB,
                DataLoc
            });

            CWChart.DataSource.AddRange(new List<ChartData1D>()
            {
                Freq1,
                Freq2,
                DataCW1,
                DataCW2
            });

            CWPoints.ItemChanged += UpdatePointPanel;
        }

        private void UpdatePointPanel(object sender, RoutedEventArgs e)
        {
            UpdateCWPointsDisplay();
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            UpdateChartAndDataFlow(true);
        }

        private void Close(object sender, RoutedEventArgs e)
        {
            Hide();
        }

        /// <summary>
        /// 刷新CW点显示
        /// </summary>
        private void UpdateCWPointsDisplay()
        {
            DataViewer.ClearItems();
            for (int i = 0; i < CWPoints.Count; i++)
            {
                DataViewer.AddItem(CWPoints[i],
                    Math.Round(CWPoints[i].MoverLoc, 4),
                    Math.Round(CWPoints[i].CW1, 4),
                    Math.Round(CWPoints[i].CW2, 4),
                    Math.Round(CWPoints[i].Bv, 4),
                    Math.Round(CWPoints[i].Bp, 4),
                    Math.Round(CWPoints[i].B, 4));
            }
        }

        /// <summary>
        /// 清除数据
        /// </summary>
        public void ClearData()
        {
            DataBp.Data.Clear();
            DataBp.Data.Clear();
            DataB.Data.Clear();
            DataCW1.Data.Clear();
            DataCW2.Data.Clear();
            Freq1.Data.Clear();
            Freq2.Data.Clear();
            DataLoc.Data.Clear();

            DataViewer.ClearItems();
            Chart.UpdateChartAndDataFlow(true);
            CWChart.UpdateChartAndDataFlow(true);
        }

        /// <summary>
        /// 数据点选中事件
        /// </summary>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        private void SelectPoint(int arg1, object arg2)
        {
            UpdateChartAndDataFlow(true);
        }

        /// <summary>
        /// 刷新CW图表
        /// </summary>
        public void UpdateChartAndDataFlow(bool isAutoScale)
        {
            DataLoc.Data = CWPoints.Select((x) => x.MoverLoc).ToList();
            DataBp.Data = CWPoints.Select((x) => x.Bp).ToList();
            DataBv.Data = CWPoints.Select((x) => x.Bv).ToList();
            DataB.Data = CWPoints.Select((x) => x.B).ToList();

            Chart.UpdateChartAndDataFlow(isAutoScale);

            if (DataViewer.GetSelectedTag() == null) return;
            CWPointObject point = DataViewer.GetSelectedTag() as CWPointObject;

            DataCW1.Data = point.CW1Values;
            DataCW2.Data = point.CW2Values;
            Freq1.Data = point.CW1Freqs;
            Freq2.Data = point.CW2Freqs;

            CWChart.UpdateChartAndDataFlow(isAutoScale);
        }

        /// <summary>
        /// 用二次函数拟合
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void FitWithCurve(object sender, RoutedEventArgs e)
        {
            if (IsAngleWindow)
            {
            }
            else
            {

            }
        }

        /// <summary>
        /// 将绝对值数据转换成三角函数
        /// </summary>
        public List<double> ConvertAbsDataToSin(List<double> y)
        {
            List<double> ny = new List<double>();
            if (y.Count < 3)
            {
                return y;
            }

            List<int> revs = new List<int>();

            for (int i = 1; i < y.Count - 1; i++)
            {
                //找拐点
                if (y[i] <= y[i + 1] && y[i] <= y[i - 1])
                {
                    revs.Add(i);
                }
            }

            revs.Sort();

            for (int i = 0; i < y.Count; i++)
            {
                int sgn = 1;
                foreach (var item in revs)
                {
                    if (i > item) sgn *= -1;
                }
                ny.Add(sgn * y[i]);
            }

            return ny;
        }

        /// <summary>
        /// 切换视图
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void ChangeView(object sender, RoutedEventArgs e)
        {
            CWChart.Visibility = Visibility.Collapsed;
            Chart.Visibility = Visibility.Collapsed;
            if (sender == MagnetBtn)
            {
                MagnetBtn.KeepPressed = true;
                CWBtn.KeepPressed = false;
                Chart.Visibility = Visibility.Visible;
            }
            if (sender == CWBtn)
            {
                MagnetBtn.KeepPressed = false;
                CWBtn.KeepPressed = true;
                CWChart.Visibility = Visibility.Visible;
            }
        }
    }
}
