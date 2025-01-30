﻿using Controls;
using HardWares.纳米位移台;
using ODMR_Lab.Windows;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Shapes;
using ComboBox = Controls.ComboBox;

namespace ODMR_Lab.位移台部分
{
    /// <summary>
    /// MoverTestWindow.xaml 的交互逻辑
    /// </summary>
    public partial class MoverTestWindow : Window
    {
        PartTypes part;

        public MoverTestWindow(PartTypes part, string title)
        {
            InitializeComponent();
            this.part = part;
            Title.Content = title;
            List<NanoStageInfo> movers = DeviceDispatcher.TryGetMoverDevice(part);
            //占用设备
            DeviceDispatcher.UseDevices(movers.Select(x => x as InfoBase).ToList());
            foreach (var item in movers)
            {
                MoverLists.Children.Add(CreateMoverGridBar(item));
            }
            CreateListenThread();
        }


        private void Close(object sender, RoutedEventArgs e)
        {
            DisposeThread();
            List<NanoStageInfo> movers = DeviceDispatcher.TryGetMoverDevice(part);
            //解除占用
            DeviceDispatcher.EndUseDevices(movers.Select(x => x as InfoBase).ToList());
            Close();
        }

        private void Window_MouseMove(object sender, MouseEventArgs e)
        {
            if (e.LeftButton == MouseButtonState.Pressed && e.GetPosition(this).Y < 30)
            {
                DragMove();
            }
        }


        #region 控制器信息条

        /// <summary>
        /// 创建位移台标签
        /// </summary>
        /// <returns></returns>
        private Grid CreateMoverGridBar(NanoStageInfo info)
        {
            Grid grid = new Grid();
            grid.Tag = info;
            grid.Height = 50;
            grid.ColumnDefinitions.Add(new ColumnDefinition());
            grid.ColumnDefinitions.Add(new ColumnDefinition() { Width = new GridLength(60) });
            grid.ColumnDefinitions.Add(new ColumnDefinition() { Width = new GridLength(90) });
            grid.ColumnDefinitions.Add(new ColumnDefinition() { Width = new GridLength(90) });
            grid.ColumnDefinitions.Add(new ColumnDefinition() { Width = new GridLength(60) });
            grid.ColumnDefinitions.Add(new ColumnDefinition() { Width = new GridLength(40) });
            grid.ColumnDefinitions.Add(new ColumnDefinition() { Width = new GridLength(40) });

            Label l = CreateLabelStyle();
            l.Content = info.Parent.Device.ProductName;
            grid.Children.Add(l);
            Grid.SetColumn(l, 0);

            l = CreateLabelStyle();
            l.Content = info.Device.AxisName;
            grid.Children.Add(l);
            Grid.SetColumn(l, 1);

            ComboBox box = CreateComboxStyle(info);
            grid.Children.Add(box);
            Grid.SetColumn(box, 2);

            FontChangeText target = new FontChangeText();
            target.Margin = new Thickness(5);
            target.InnerTextBox.IsReadOnly = true;
            target.Tag = info;
            TextTemplate.CloneStyleTo(target);
            grid.Children.Add(target);
            Grid.SetColumn(target, 3);

            FontChangeText text = new FontChangeText();
            text.InnerTextBox.Text = "0.01";
            text.Margin = new Thickness(5);
            TextTemplate.CloneStyleTo(text);
            grid.Children.Add(text);
            Grid.SetColumn(text, 4);

            DecoratedButton btn = new DecoratedButton() { Text = "<" };
            BtnTemplate.CloneStyleTo(btn);
            btn.Tag = info;
            btn.Click += MoveLeft;
            grid.Children.Add(btn);
            btn.Margin = new Thickness(5);
            btn.CornerRadius = new CornerRadius(5);
            Grid.SetColumn(btn, 5);

            btn = new DecoratedButton() { Text = ">" };
            BtnTemplate.CloneStyleTo(btn);
            btn.Tag = info;
            btn.Click += MoveRight;
            grid.Children.Add(btn);
            btn.Margin = new Thickness(5);
            btn.CornerRadius = new CornerRadius(5);
            Grid.SetColumn(btn, 6);

            return grid;
        }


        private Label CreateLabelStyle()
        {
            Label l = new Label();
            l.Foreground = Brushes.White;
            l.FontSize = 10;
            l.HorizontalAlignment = HorizontalAlignment.Stretch;
            l.VerticalAlignment = VerticalAlignment.Stretch;
            l.HorizontalContentAlignment = HorizontalAlignment.Center;
            l.VerticalContentAlignment = VerticalAlignment.Center;
            return l;
        }

        private ComboBox CreateComboxStyle(NanoStageInfo info)
        {
            ComboBox box = new ComboBox();
            BtnTemplate.CloneStyleTo(box);
            box.FontSize = 12;
            box.TextAreaRatio = BtnTemplate.TextAreaRatio;
            box.IconSource = BtnTemplate.IconSource;
            box.ImagePlace = BtnTemplate.ImagePlace;
            box.PanelWidth = 200;
            box.IconMargin = BtnTemplate.IconMargin;
            box.TemplateButton = BtnTemplate;
            foreach (var item in Enum.GetNames(typeof(MoverTypes)))
            {
                box.Items.Add(new DecoratedButton() { Text = item });
            }
            box.Select(Enum.GetName(info.MoverType.GetType(), info.MoverType));
            box.Margin = new Thickness(5);
            return box;
        }
        #endregion


        #region 位置监听线程

        Thread listenthread = null;
        bool IsthreadEnd = false;
        private void CreateListenThread()
        {
            listenthread = new Thread(() =>
            {
                while (!IsthreadEnd)
                {
                    try
                    {
                        List<NanoStageInfo> movers = DeviceDispatcher.TryGetMoverDevice(part);
                        foreach (var item in movers)
                        {
                            double pos = item.Device.Position;
                            Dispatcher.Invoke(() =>
                            {
                                FontChangeText text = GetTextbox(item);
                                if (text != null)
                                {
                                    text.InnerTextBox.Text = pos.ToString();
                                }
                            });
                        }

                        Thread.Sleep(100);
                    }
                    catch (Exception ex) { }
                }
            });
            listenthread.Start();
        }

        private FontChangeText GetTextbox(NanoStageInfo mover)
        {
            foreach (var item in MoverLists.Children)
            {
                if ((item as Grid).Tag as NanoStageInfo == mover)
                {
                    return ((item as Grid).Children[4] as FontChangeText);
                }
            }
            return null;
        }

        private void DisposeThread()
        {
            IsthreadEnd = true;
            while (listenthread.ThreadState == ThreadState.Running)
            {
                Thread.Sleep(50);
            }
        }
        #endregion

        #region 移动位移台

        Thread MoveThread = null;
        bool IsMoveStopped = false;

        private void MoveLeft(object sender, RoutedEventArgs e)
        {
            bool result = double.TryParse((((sender as DecoratedButton).Parent as Grid).Children[5] as FontChangeText).InnerTextBox.Text, out double step);
            NanoStageInfo info = (sender as DecoratedButton).Tag as NanoStageInfo;
            if (result == false) return;
            StageBase stage = info.Device;
            if (MoveThread == null)
            {
                MoveThread = new Thread(() =>
                {
                    IsMoveStopped = false;
                    stage.MoveStepAndWait(-step, 5000);
                    IsMoveStopped = true;
                });
                MoveThread.Start();
                return;
            }
            if (IsMoveStopped == false)
            {
                return;
            }
            MoveThread = new Thread(() =>
            {
                IsMoveStopped = false;
                stage.MoveStepAndWait(-step, 5000);
                IsMoveStopped = true;
            });
            MoveThread.Start();
            return;
        }

        private void MoveRight(object sender, RoutedEventArgs e)
        {
            bool result = double.TryParse((((sender as DecoratedButton).Parent as Grid).Children[5] as FontChangeText).InnerTextBox.Text, out double step);
            NanoStageInfo info = (sender as DecoratedButton).Tag as NanoStageInfo;
            if (result == false) return;
            StageBase stage = info.Device;
            if (MoveThread == null)
            {
                MoveThread = new Thread(() =>
                {
                    IsMoveStopped = false;
                    stage.MoveStepAndWait(step, 5000);
                    IsMoveStopped = true;
                });
                MoveThread.Start();
                return;
            }
            if (IsMoveStopped == false)
            {
                return;
            }
            MoveThread = new Thread(() =>
            {
                IsMoveStopped = false;
                stage.MoveStepAndWait(step, 5000);
                IsMoveStopped = true;
            });
            MoveThread.Start();
            return;
        }

        #endregion

        /// <summary>
        /// 应用修改
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        private void Apply(object sender, RoutedEventArgs e)
        {
            List<MoverTypes> types = new List<MoverTypes>();
            List<NanoStageInfo> infos = new List<NanoStageInfo>();
            for (int i = 1; i < MoverLists.Children.Count; i++)
            {
                Grid g = MoverLists.Children[i] as Grid;
                infos.Add(g.Tag as NanoStageInfo);
                types.Add((MoverTypes)Enum.Parse(typeof(MoverTypes), (g.Children[2] as ComboBox).Text));
            }
            for (int i = 0; i < types.Count; i++)
            {
                MoverTypes type1 = types[i];
                for (int j = 0; j < types.Count; j++)
                {
                    MoverTypes type2 = types[j];
                    if (i != j && type1 == type2 && type1 != MoverTypes.None)
                    {
                        MessageWindow.ShowTipWindow("轴标签不能相同：" + Enum.GetName(type1.GetType(), type1), this);
                        return;
                    }
                }
            }
            for (int i = 0; i < types.Count; i++)
            {
                infos[i].MoverType = types[i];
            }
            MessageWindow.ShowTipWindow("轴参数设置成功", this);
            return;
        }
    }
}
